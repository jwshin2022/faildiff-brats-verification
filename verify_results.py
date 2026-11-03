"""Independent numerical review of reduced statistics; no image/model readers.

Standalone numerical verification code. No original model or imaging software is bundled.
Run: python reproduce_minimum.py --data minimum_data.json --out result.json
Smoke: python reproduce_minimum.py --synthetic --out smoke.json
"""
import argparse
import collections
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import resource
import sys
import time

REGIONS = ((1, 2, 3), (1, 3), (3,), (1,))
NAMES = ('WT', 'TC', 'ET', 'class1')
ARMS = ('A', 'B', 'C', 'point_B')


def dice(p, g, t):
    assert all(isinstance(a, int) and a >= 0 for a in (p, g, t))
    assert t <= min(p, g)
    return Fraction(2*t, p+g) if p+g else Fraction(1)


def scalar_condition(group, lm, scale_type, baseline_type):
    d = group['estimated_baseline_dice']
    lo_d = [max(0., a-.05*lm) for a in d]
    hi_d = [min(1., a+.05*lm) for a in d]
    center = [(a+b)/2 for a, b in zip(lo_d, hi_d)]
    radius = [(b-a)/2 for a, b in zip(lo_d, hi_d)]
    if baseline_type == 'NO_ABS_PENALTY':
        radius = [0.]*4
    if baseline_type == 'POINT_BASELINE':
        center, radius = d, [0.]*4
    counts = collections.Counter()
    full_by_region = [0]*4
    changed_passes = [0]*4
    unchanged_passes = [0]*4
    any_not_strict = False
    all_nonpositive = True
    for item in group['candidates']:
        m = item['mass']
        assert m > 0
        a0 = math.fsum(item['alpha'])
        assert a0 > 0 and all(a > 0 for a in item['alpha'])
        hs, losses, fulls, bind = [], [], [], []
        for r, labels in enumerate(REGIONS):
            p = math.fsum(item['alpha'][k] for k in labels)/a0
            ss = math.sqrt(p*(1-p)*(m+a0)/(m*(1+a0)))
            floor = 1/m
            s = {'ORIGINAL':max(ss, floor), 'STAT_ONLY':ss, 'FLOOR_ONLY':floor}[scale_type]
            low = m*max(0., p-lm*s)
            high = m*min(1., p+lm*s)
            full = low == 0 and high == m
            sf = p-lm*ss <= 0 and p+lm*ss >= 1
            ff = p-lm*floor <= 0 and p+lm*floor >= 1
            sign = int(item['target'] in labels)-int(item['source'] in labels)
            e = sign*m
            overlap = low if sign > 0 else (-high if sign < 0 else 0.)
            h = 2*overlap-center[r]*e
            single = h-radius[r]*abs(e)
            hs.append(h); losses.append(single); fulls.append(full); bind.append(floor >= ss)
            full_by_region[r] += full
            changed_passes[r] += single >= 0 and sign != 0
            unchanged_passes[r] += single >= 0 and sign == 0
            for name, yes in [('stat_only_full_entries',sf),('floor_only_full_entries',ff),
                              ('both_sufficient_entries',sf and ff),
                              ('stat_exclusive_entries',sf and not ff),
                              ('floor_exclusive_entries',ff and not sf)]:
                counts[name] += yes
        counts['candidate_count'] += 1
        counts['full_entries'] += sum(fulls)
        counts['all4_full_candidates'] += all(fulls)
        counts['singleton_pass_count'] += min(losses) >= 0
        counts['mass_le_lambda'] += m <= lm
        counts['floor_binding_entries'] += sum(bind)
        counts['floor_all4_binding'] += all(bind)
        # Same finite float64 coefficients, distinct scalar arithmetic implementation.
        any_not_strict |= sum(hs) >= 0
        all_nonpositive &= max(hs) <= 0
    counts['sum_h_impossibility_sufficient'] = not any_not_strict
    counts['all_h_nonpositive'] = all_nonpositive
    counts['full_by_region'] = full_by_region
    counts['region_changed_membership_singleton_pass_counts'] = changed_passes
    counts['region_unchanged_membership_singleton_pass_counts'] = unchanged_passes
    return dict(counts)


def reproduce(data):
    import numpy as np
    groups = data['groups']; n = len(groups)
    assert n > 0 and len({x['group_key'] for x in groups}) == n
    out = {'groups':n, 'verified_independent_patient_count':None, 'arms':{}, 'conditions':{}}
    candidates = beneficial = positive = harmful = saturated = proof = certified = 0
    integer_checks = 0
    all_deltas = {a:[] for a in ARMS}
    all_dice = {a:[] for a in ARMS}
    for group in groups:
        base = group['arms']['A']['counts']
        bd = [dice(*r) for r in base]
        for arm in ARMS:
            cc = group['arms'][arm]['counts']
            assert [x[1] for x in cc] == [x[1] for x in base]
            dd = [dice(*x) for x in cc]
            all_dice[arm].append([float(x) for x in dd])
            all_deltas[arm].append([float(x-y) for x,y in zip(dd,bd)])
            tr = group['arms'][arm]['transitions']
            assert tr['changed'] == tr['repaired']+tr['new_errors']+tr['wrong_to_wrong']
            integer_checks += 4
        for item in group['candidates']:
            delta = []
            for r, labels in enumerate(REGIONS):
                p,g,t = base[r]
                e = item['mass']*(int(item['target'] in labels)-int(item['source'] in labels))
                dt = item['gt_tp_increment'][r]
                dd = dice(p+e,g,t+dt)-bd[r]
                delta.append(dd)
                integer_checks += 1
            u = sum(delta[:3])/3
            h = min(delta) < 0
            candidates += 1; positive += u > 0; harmful += h
            beneficial += u > 0 and not h
        orig = scalar_condition(group,data['lambda_original'],'ORIGINAL','ORIGINAL')
        saturated += orig['full_entries']; proof += orig['sum_h_impossibility_sufficient']
        certified += orig['singleton_pass_count']
    seed, draws_n = data['bootstrap']['seed'], data['bootstrap']['draws']
    draws = np.random.Generator(np.random.PCG64(seed)).integers(0,n,size=(draws_n,n))
    for arm in ARMS:
        dd = np.asarray(all_deltas[arm],dtype=np.float64)
        dm = dd[:,:3].mean(axis=1)
        means = np.asarray(all_dice[arm],dtype=np.float64)
        out['arms'][arm] = {
            'region_mean_Dice':dict(zip(NAMES,means.mean(axis=0).tolist())),
            'region_mean_delta':dict(zip(NAMES,dd.mean(axis=0).tolist())),
            'mean_Dice':float(means[:,:3].mean(axis=1).mean()),
            'mean_delta':float(dm.mean()),
            'mean_delta_CI95':np.quantile(dm[draws].mean(axis=1),[.025,.975],method='linear').tolist(),
            'region_delta_CI95':{r:np.quantile(dd[:,k][draws].mean(axis=1),[.025,.975],method='linear').tolist() for k,r in enumerate(NAMES)},
            'acted':sum(g['arms'][arm]['transitions']['changed']>0 for g in groups),
            'any_region_harm':int((dd.min(axis=1)<0).sum()),
            'improved':int((dm>0).sum()), 'worsened':int((dm<0).sum()),
            'worst_mean_delta':float(dm.min()),
            'transitions':{k:sum(g['arms'][arm]['transitions'][k] for g in groups) for k in ('changed','repaired','new_errors','wrong_to_wrong')}}
    out['mechanism'] = {'candidates':candidates,'safe_beneficial_singletons':beneficial,
                        'positive_mean_singletons':positive,'any_region_harm_singletons':harmful,
                        'saturated_intervals':saturated,'region_intervals':candidates*4,
                        'joint_exclusion_groups':proof,'certificate_singleton_pass':certified}
    for lm_name,lm in [('ORIGINAL',data['lambda_original']),('CELL_COMPONENT_REFERENCE',data['lambda_cell_reference'])]:
        for sc in ('ORIGINAL','STAT_ONLY','FLOOR_ONLY'):
            for bs in ('ORIGINAL','NO_ABS_PENALTY','POINT_BASELINE'):
                rows = [scalar_condition(g,lm,sc,bs) for g in groups]
                fields = set(rows[0])-{'full_by_region','region_changed_membership_singleton_pass_counts','region_unchanged_membership_singleton_pass_counts'}
                acc = {k:sum(x[k] for x in rows) for k in fields}
                for k in ('full_by_region','region_changed_membership_singleton_pass_counts','region_unchanged_membership_singleton_pass_counts'):
                    acc[k] = [sum(x[k][r] for x in rows) for r in range(4)]
                out['conditions']['/'.join((lm_name,sc,bs))] = acc
    out['integer_count_checks'] = integer_checks
    out['bootstrap'] = dict(data['bootstrap'],generator='PCG64',quantile='linear',
                             interpretation='exploratory source/image-group resampling; patient independence unverified')
    out['environment'] = {'python':platform.python_version(),'numpy':np.__version__,'platform':platform.platform()}
    return out


def synthetic():
    # Fully artificial statistics, not any patient or source dataset.
    assert dice(0,0,0)==1 and dice(0,2,0)==0 and dice(2,2,1)==Fraction(1,2)
    baseline = [[2,3,2],[1,2,1],[1,1,1],[0,1,0]]
    arm = {'counts':baseline,'transitions':{'changed':0,'repaired':0,'new_errors':0,'wrong_to_wrong':0}}
    row = {'group_key':'synthetic-01','arms':{a:arm for a in ARMS},
           'estimated_baseline_dice':[.5]*4,
           'candidates':[{'mass':1,'source':2,'target':1,'gt_tp_increment':[0,1,0,1],
                          'alpha':[1.,1.,1.,1.]}]}
    data={'groups':[row], 'lambda_original':17.399819144085278,
          'lambda_cell_reference':6.977041320835106,'bootstrap':{'seed':20261008,'draws':100}}
    result=reproduce(data)
    assert result['mechanism']['safe_beneficial_singletons']==1
    assert result['mechanism']['saturated_intervals']==4
    assert result['mechanism']['joint_exclusion_groups']==1
    assert len(result['conditions'])==18
    assert result['arms']['A']['mean_delta_CI95']==[0.,0.]
    return {'status':'PASS','tests':['both_empty','one_empty','rational_Dice','nested_region_increment',
             'safe_beneficial_candidate','full_saturation','joint_sufficient_exclusion','18_conditions','zero_bootstrap'],
            'patient_data':False,'numpy_version':result['environment']['numpy']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--data',type=Path);p.add_argument('--out',type=Path,required=True);p.add_argument('--synthetic',action='store_true')
    a=p.parse_args();started=time.process_time();wall=time.monotonic()
    result=synthetic() if a.synthetic else reproduce(json.loads(a.data.read_text()))
    result['execution']={'cpu_seconds':time.process_time()-started,'wall_seconds':time.monotonic()-wall,
                         'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                         'raw_image_reads':0,'GT_array_reads':0,'GPU_used':False,'production_analysis_imports':False}
    with a.out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps({'status':'PASS','mode':'synthetic' if a.synthetic else 'reduced-derived',
                      'groups':result.get('groups'),'cpu_seconds':result['execution']['cpu_seconds']}))


if __name__=='__main__':
    try:main()
    except Exception as e:
        print(json.dumps({'status':'FAIL_CLOSED','type':type(e).__name__}));sys.exit(2)
