"""Bounded reviewer-requested analyses from count summaries; no model/image loader.

All three confidence thresholds are reported. No threshold is fitted to outcomes.
"""
import argparse, json, hashlib, time, platform
from pathlib import Path
from fractions import Fraction
import numpy as np

REGIONS=((1,2,3),(1,3),(3,),(1,))
def dice(c):
    p,g,t=map(int,c)
    assert p>=0 and g>=0 and 0<=t<=min(p,g)
    return Fraction(2*t,p+g) if p+g else Fraction(1)
def interval(x,stat,seed=20261008):
    x=np.asarray(x,float);rng=np.random.Generator(np.random.PCG64(seed))
    ix=rng.integers(0,len(x),size=(20000,len(x)))
    z=np.mean(x[ix],axis=1) if stat=='mean' else np.median(x[ix],axis=1)
    return np.quantile(z,[.025,.975],method='linear').tolist()
def describe(x):
    a=np.asarray(x,float)
    return dict(n=len(a),mean=float(a.mean()),median=float(np.median(a)),
                q1=float(np.quantile(a,.25)),q3=float(np.quantile(a,.75)),
                minimum=float(a.min()),maximum=float(a.max()))
def evaluate(group,selected):
    base=group['arms']['A']['counts'];counts=[list(z) for z in base]
    tr=dict(changed=0,repaired=0,new_errors=0,wrong_to_wrong=0)
    for c in selected:
        m=c['mass'];tr['changed']+=m
        for k in ('repaired','new_errors','wrong_to_wrong'):tr[k]+=c['transitions'][k]
        for r,labels in enumerate(REGIONS):
            counts[r][0]+=m*(int(c['target'] in labels)-int(c['source'] in labels))
            counts[r][2]+=c['gt_tp_increment'][r]
    assert tr['changed']==tr['repaired']+tr['new_errors']+tr['wrong_to_wrong']
    delta=[dice(c)-dice(b) for c,b in zip(counts,base)]
    return dict(delta=[float(z) for z in delta],utility=float(sum(delta[:3])/3),
                any_harm=min(delta)<0,selected=len(selected),transitions=tr)
def run(d):
    groups=d['groups'];assert len(groups)==93
    out={'groups':len(groups),'bootstrap':{'draws':20000,'generator':'PCG64',
         'full_sample_seed':20261008,'sensitivity_seed':20261103},'arms':{},'confidence_baselines':{},'checks':{}}
    utilities={};paired=[]
    for arm in ('B','point_B'):
        rows=[];ratios=[];inf=zerozero=0
        for g in groups:
            delta=[float(dice(x)-dice(b)) for x,b in zip(g['arms'][arm]['counts'],g['arms']['A']['counts'])]
            rows.append(sum(delta[:3])/3);t=g['arms'][arm]['transitions'];r,n=t['repaired'],t['new_errors']
            if r:ratios.append(n/r)
            elif n:inf+=1
            else:zerozero+=1
        utilities[arm]=np.asarray(rows)
        out['arms'][arm]=dict(utility=describe(rows),mean_CI95=interval(rows,'mean'),
          median_CI95=interval(rows,'median'),improved=sum(x>0 for x in rows),worsened=sum(x<0 for x in rows),
          ratio_finite=describe(ratios) if ratios else None,ratio_infinite=inf,ratio_undefined=zerozero,
          groups_new_exceeds_repaired=sum(g['arms'][arm]['transitions']['new_errors']>g['arms'][arm]['transitions']['repaired'] for g in groups),
          groups_equal=sum(g['arms'][arm]['transitions']['new_errors']==g['arms'][arm]['transitions']['repaired'] for g in groups))
    worst=np.argsort(utilities['B'],kind='stable')
    out['sensitivity']={}
    for k in (1,3):
        keep=np.ones(len(groups),bool);keep[worst[:k]]=False
        out['sensitivity'][str(k)]={'selection':'post hoc worst B utility; same exclusions for both arms','remaining':int(keep.sum()),
          'excluded_B_utility':utilities['B'][worst[:k]].tolist(),
          'excluded_sum_share_of_full_sum':float(utilities['B'][worst[:k]].sum()/utilities['B'].sum()),
          'arms':{a:dict(**describe(v[keep]),mean_CI95=interval(v[keep],'mean',20261103)) for a,v in utilities.items()}}
    allc=[c for g in groups for c in g['candidates']];m=np.array([c['mass'] for c in allc])
    q=np.array([np.asarray(c['alpha'])/sum(c['alpha']) for c in allc]);b=np.array([c['source'] for c in allc]);t=np.array([c['target'] for c in allc])
    margin=q[np.arange(len(q)),t]-q[np.arange(len(q)),b];score=m*margin
    out['candidate_generation']=dict(mass=describe(m),one_voxel=int((m==1).sum()),mass_le_lambda=int((m<=d['lambda_original']).sum()),
      margin=describe(margin),score=describe(score),positive_margin=int((margin>0).sum()),zero_margin=int((margin==0).sum()),negative_margin=int((margin<0).sum()),
      target_probability=describe(q[np.arange(len(q)),t]))
    # Exact aggregation of disjoint cells validated against the stored completed B arm.
    checks=0
    for g in groups:
        row=evaluate(g,g['candidates'])
        expected=[float(dice(c)-dice(a)) for c,a in zip(g['arms']['B']['counts'],g['arms']['A']['counts'])]
        assert row['delta']==expected and row['transitions']==g['arms']['B']['transitions'];checks+=1
    out['checks']['all_candidates_reproduce_B']=checks
    for threshold in (.5,.7,.9):
        rows=[]
        for g in groups:
            picks=[]
            for c in g['candidates']:
                z=np.asarray(c['alpha'])/sum(c['alpha'])
                if z[c['target']]>=threshold and z[c['target']]>z[c['source']]:picks.append(c)
            rows.append(evaluate(g,picks))
        u=[x['utility'] for x in rows]
        out['confidence_baselines'][str(threshold)]=dict(utility=describe(u),mean_CI95=interval(u,'mean'),
          region_mean_delta=np.mean([x['delta'] for x in rows],axis=0).tolist(),
          acted=sum(x['selected']>0 for x in rows),harm=sum(x['any_harm'] for x in rows),
          improved=sum(v>0 for v in u),worsened=sum(v<0 for v in u),
          selected=sum(x['selected'] for x in rows),
          transitions={k:sum(x['transitions'][k] for x in rows) for k in rows[0]['transitions']})
    out['design_rules']={'lambda':d['lambda_original'],'minimum_integer_mass_potentially_not_floor_saturated':18,
      'minimum_integer_mass_guaranteed_not_floor_saturated':35,
      'quantile_n32':32,'minimum_n_nonmaximum_95pct':39,
      'joint_four_bonferroni_example_min_n_nonmaximum':159,
      'note':'Bonferroni example is not the original maximum-score calibration.'}
    # Plot points: no identity columns; pair retained only within same anonymous row.
    for i,g in enumerate(groups):
        item={a:dict(utility=float(utilities[a][i]),**g['arms'][a]['transitions']) for a in utilities};paired.append(item)
    return out,paired,worst[:3].tolist()
def main():
    p=argparse.ArgumentParser();p.add_argument('--data',required=True);p.add_argument('--out',required=True);a=p.parse_args()
    root=Path(a.out);root.mkdir(parents=True,exist_ok=False);start=time.time();raw=Path(a.data).read_bytes();d=json.loads(raw)
    out,paired,worst=run(d);out['environment']={'python':platform.python_version(),'numpy':np.__version__,'seconds':time.time()-start}
    (root/'REVISION_ANALYSES.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    (root/'GROUP_FIGURE_INPUT_PRIVATE.json').write_text(json.dumps(paired,indent=2)+'\n')
    (root/'PRIVATE_PROVENANCE.json').write_text(json.dumps({'input_sha256':hashlib.sha256(raw).hexdigest(),'worst_indices':worst,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
