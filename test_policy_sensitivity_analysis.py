import itertools,json
from fractions import Fraction
from policy_sensitivity_analysis import composite,REGIONS
def counts(y,p):return [[sum(v in ls for v in p),sum(v in ls for v in y),sum(a in ls and b in ls for a,b in zip(p,y))] for ls in REGIONS]
def ratio(c):return Fraction(2*c[2],c[0]+c[1]) if c[0]+c[1] else Fraction(1)
n=0
for a in itertools.product(range(4),repeat=2):
 for b in itertools.product(range(4),repeat=2):
  for y in itertools.product(range(4),repeat=2):
   base=counts(y,a);cs=[]
   for j in range(2):
    if b[j]==a[j]:continue
    single=list(a);single[j]=b[j];cc=counts(y,single)
    cs.append((j,{'mass':1,'source':a[j],'target':b[j],'gt_tp_increment':[c[2]-v[2] for c,v in zip(cc,base)],
      'transitions':{'repaired':int(a[j]!=y[j] and b[j]==y[j]),'new_errors':int(a[j]==y[j] and b[j]!=y[j]),'wrong_to_wrong':int(a[j]!=y[j] and b[j]!=y[j])}}))
   for mask in itertools.product([0,1],repeat=len(cs)):
    pred=list(a);selected=[]
    for (j,c),take in zip(cs,mask):
     if take:pred[j]=b[j];selected.append(c)
    expected=counts(y,pred);d=[ratio(c)-ratio(v) for c,v in zip(expected,base)]
    z=composite({'arms':{'A':{'counts':base}}},selected)
    assert z['counts']==expected and z['delta_pp']==[float(x*100) for x in d]
    assert z['utility_pp']==float(sum(d[:3])*100/3) and z['harm']==any(x<0 for x in d)
    assert z['transitions']=={'changed':sum(x!=v for x,v in zip(pred,a)),
       'repaired':sum(x!=v and x==t for x,v,t in zip(pred,a,y)),
       'new_errors':sum(x!=v and v==t for x,v,t in zip(pred,a,y)),
       'wrong_to_wrong':sum(x!=v and x!=t and v!=t for x,v,t in zip(pred,a,y))}
    n+=1
print(json.dumps({'status':'PASS','direct_array_vs_integer_composite_checks':n,'base_label_combinations':4096,'empty_regions_included':True}))
