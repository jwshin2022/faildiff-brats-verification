"""Independent exhaustive array checks for disjoint count aggregation."""
import itertools,json
from fractions import Fraction
from revision_analysis import evaluate,REGIONS
def counts(y,p):return [[sum(v in ls for v in p),sum(v in ls for v in y),sum(a in ls and b in ls for a,b in zip(p,y))] for ls in REGIONS]
def ratio(c):return Fraction(2*c[2],c[0]+c[1]) if c[0]+c[1] else Fraction(1)
a=[0,1,2,3];target=[2,3,1,0];tests=0
for y in itertools.product(range(4),repeat=4):
 base=counts(y,a);candidates=[]
 for j in range(4):
  p=a.copy();p[j]=target[j];cc=counts(y,p)
  candidates.append(dict(mass=1,source=a[j],target=target[j],gt_tp_increment=[x[2]-b[2] for x,b in zip(cc,base)],
   transitions={'repaired':int(a[j]!=y[j] and p[j]==y[j]),'new_errors':int(a[j]==y[j] and p[j]!=y[j]),'wrong_to_wrong':int(a[j]!=y[j] and p[j]!=y[j])}))
 g={'arms':{'A':{'counts':base}}}
 for mask in itertools.product([0,1],repeat=4):
  p=[target[j] if mask[j] else a[j] for j in range(4)];cc=counts(y,p);delta=[ratio(x)-ratio(b) for x,b in zip(cc,base)]
  z=evaluate(g,[c for c,k in zip(candidates,mask) if k])
  assert z['delta']==[float(x) for x in delta] and z['utility']==float(sum(delta[:3])/3)
  assert z['any_harm']==(min(delta)<0)
  assert z['transitions']=={'changed':sum(mask),'repaired':sum(v!=w and w==t for v,w,t in zip(a,p,y)),
        'new_errors':sum(v!=w and v==t for v,w,t in zip(a,p,y)),
        'wrong_to_wrong':sum(v!=w and v!=t and w!=t for v,w,t in zip(a,p,y))};tests+=1
print(json.dumps({'status':'PASS','synthetic_array_vs_count_checks':tests,'empty_region_cases_included':True,'study_data_checked_by_this_test':False}))
