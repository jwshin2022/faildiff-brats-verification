"""Policy sensitivity analyses from fixed integer summaries; no image/model access."""
from pathlib import Path
from fractions import Fraction
import argparse, json, hashlib, platform, time
import numpy as np

REGIONS=((1,2,3),(1,3),(3,),(1,))
NAMES=('WT','TC','ET','NCR')
def dice(c):
    p,g,t=map(int,c)
    assert min(p,g,t)>=0 and t<=min(p,g)
    return Fraction(2*t,p+g) if p+g else Fraction(1)

def composite(group,selected):
    base=group['arms']['A']['counts']
    counts=[list(map(int,v)) for v in base]
    tr={k:0 for k in ('changed','repaired','new_errors','wrong_to_wrong')}
    for c in selected:
        m=int(c['mass']);tr['changed']+=m
        for k in ('repaired','new_errors','wrong_to_wrong'): tr[k]+=int(c['transitions'][k])
        for r,labels in enumerate(REGIONS):
            counts[r][0]+=m*(int(c['target'] in labels)-int(c['source'] in labels))
            counts[r][2]+=int(c['gt_tp_increment'][r])
    delta=[dice(a)-dice(b) for a,b in zip(counts,base)]
    assert tr['changed']==sum(tr[k] for k in ('repaired','new_errors','wrong_to_wrong'))
    return {'counts':counts,'baseline_counts':base,'delta_pp':[float(x*100) for x in delta],
            'utility_pp':float(sum(delta[:3])*100/3),'harm':any(x<0 for x in delta),
            'selected':len(selected),'transitions':tr}

def bootstrap(x,median=False,seed=20261008):
    x=np.asarray(x,float);rng=np.random.Generator(np.random.PCG64(seed))
    ix=rng.integers(0,len(x),size=(20000,len(x)))
    z=np.median(x[ix],axis=1) if median else np.mean(x[ix],axis=1)
    return np.quantile(z,[.025,.975],method='linear').tolist()

def summary(rows):
    u=np.array([r['utility_pp'] for r in rows]);a=np.array([r['selected']>0 for r in rows])
    regions=np.array([r['delta_pp'] for r in rows]);worst=np.argsort(u,kind='stable')
    tr={k:sum(r['transitions'][k] for r in rows) for k in rows[0]['transitions']}
    out={'n':len(rows),'selected':sum(r['selected'] for r in rows),'acted':int(a.sum()),
         'harm':sum(r['harm'] for r in rows),'improved':int((u>0).sum()),'worsened':int((u<0).sum()),
         'unchanged':int((u==0).sum()),'mean_pp':float(u.mean()),'median_pp':float(np.median(u)),
         'iqr_pp':np.quantile(u,[.25,.75]).tolist(),'mean_CI95_pp':bootstrap(u),
         'median_CI95_pp':bootstrap(u,True),'minimum_pp':float(u.min()),'maximum_pp':float(u.max()),
         'edited_mean_pp':float(u[a].mean()),'edited_median_pp':float(np.median(u[a])),
         'transitions':tr,'region_mean_pp':dict(zip(NAMES,regions.mean(0).tolist())),
         'region_median_pp':dict(zip(NAMES,np.median(regions,axis=0).tolist())),
         'region_worsened_groups':dict(zip(NAMES,(regions<0).sum(0).tolist())),
         'region_mean_CI95_pp':{n:bootstrap(regions[:,i]) for i,n in enumerate(NAMES)},
         'sensitivity':{},'worst_three':[]}
    for k in (1,3):
        keep=np.ones(len(rows),bool);keep[worst[:k]]=False
        out['sensitivity'][str(k)]={'remaining':int(keep.sum()),'mean_pp':float(u[keep].mean()),
            'median_pp':float(np.median(u[keep])),'mean_CI95_pp':bootstrap(u[keep],seed=20261104),
            'excluded_sum_pp':float(u[worst[:k]].sum()),
            'excluded_share_of_signed_sum':float(u[worst[:k]].sum()/u.sum()) if u.sum()!=0 else None}
    for idx in worst[:3]:
        r=rows[int(idx)]
        out['worst_three'].append({k:r[k] for k in ('utility_pp','delta_pp','selected','transitions','counts','baseline_counts')})
    return out

def run(data):
    groups=data['groups'];assert len(groups)==93
    results={};private={};verified=0;positive=negative=zero=0
    for g in groups:
        z=composite(g,g['candidates']);b=g['arms']['B']
        assert z['counts']==b['counts'] and z['transitions']==b['transitions'];verified+=1
        for c in g['candidates']:
            diff=c['alpha'][c['target']]-c['alpha'][c['source']]
            positive+=int(diff>0);negative+=int(diff<0);zero+=int(diff==0)
    assert (positive,negative,zero)==(19927,3881,0)
    for name,tau in [('confidence_0.5',.5),('confidence_0.7',.7),('confidence_0.9',.9),('B_positive',None)]:
        rows=[]
        for g in groups:
            selected=[]
            for c in g['candidates']:
                q=np.asarray(c['alpha'],float);q=q/q.sum()
                if q[c['target']]>q[c['source']] and (tau is None or q[c['target']]>=tau):selected.append(c)
            rows.append(composite(g,selected))
        private[name]=rows;results[name]=summary(rows)
    assert results['B_positive']['selected']==19927
    return {'groups':93,'regions':NAMES,'score_counts':{'positive':positive,'negative':negative,'zero':zero},
            'B_complete_integer_reproduction_groups':verified,'bootstrap':{'draws':20000,'generator':'PCG64',
            'full_seed':20261008,'post_hoc_tail_seed':20261104},'policies':results},private

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--root',required=True);a=ap.parse_args()
    r=Path(a.root);r.mkdir(parents=True,exist_ok=True)
    (r/'analysis').mkdir(exist_ok=True);(r/'private').mkdir(exist_ok=True,mode=0o700)
    start=time.process_time();raw=Path(a.data).read_bytes();d=json.loads(raw)
    public,private=run(d)
    for policy in public['policies'].values():
        for row in policy['worst_three']:
            row.pop('counts');row.pop('baseline_counts')
    public['execution']={'CPU_seconds':time.process_time()-start,'Python':platform.python_version(),'NumPy':np.__version__,
        'MRI_GT_array_reads':0,'checkpoint_reads':0,'GPU_used':False}
    for name,obj in [('analysis/POLICY_SENSITIVITY_RESULTS.json',public),('private/POLICY_GROUP_ROWS.json',private),
                     ('private/ANALYSIS_BINDING.json',{'input_sha256':hashlib.sha256(raw).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})]:
        with (r/name).open('x') as f:json.dump(obj,f,indent=2,allow_nan=False)
    print(json.dumps({k:{z:v[z] for z in ('selected','acted','mean_pp','median_pp','edited_mean_pp','region_mean_pp','transitions','sensitivity')} for k,v in public['policies'].items()},indent=2))
if __name__=='__main__':main()
