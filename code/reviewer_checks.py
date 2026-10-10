"""Reviewer checks written independently of the submission's diagnostic.

Exact algebra uses Fraction and the field Q(sqrt(3)). Optional quadrature uses
the Rayleigh radial stop-loss identity followed by one angular integration,
not the submission's conditional half-normal integrand. Numerical results are
corroboration only: scipy's estimated errors are not certified bounds.
"""
from fractions import Fraction as F
import argparse
from pathlib import Path
import json
import math
import platform
import scipy
from scipy.integrate import quad
from scipy.special import erfc

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--submission', type=Path, default=Path.cwd() / 'submission',
                    help='Directory containing the extracted submission (default: ./submission).')
args = parser.parse_args()
ROOT = args.submission.expanduser().resolve()
data = json.loads((ROOT / 'data/counterexample.json').read_text())
weights = [F(s) for s in data['weights']]
sigma = [[F(v) for v in row] for row in data['target_covariance']]
components = [[[F(v) for v in row] for row in cov]
              for cov in data['component_covariances']]
vectors = [(F(4,9),F(0)), (F(4,3),F(4,3)), (F(4,3),F(-4,3))]
checks = []

def check(name, truth, evidence):
    if not truth:
        raise ValueError(name)
    checks.append({'name':name, 'passed':True, 'evidence':evidence})

def outer(v):
    return [[v[i]*v[j] for j in range(2)] for i in range(2)]

def det(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]

check('weights positive and sum exactly one',
      all(p>0 for p in weights) and sum(weights)==1,
      [str(sum(weights)),*[str(p) for p in weights]])
check('target positive definite', sigma[0][0]>0 and det(sigma)>0,
      {'determinant':str(det(sigma))})
for i, (v, cov) in enumerate(zip(vectors, components)):
    check(f'component {i} is exact rank-one outer product',
          outer(v)==cov and det(cov)==0 and sum(cov[j][j] for j in range(2))>0,
          {'vector':[str(x) for x in v], 'determinant':str(det(cov))})
check('weighted standard-deviation linear coefficients',
      [weights[0]*vectors[0][0],weights[1]*vectors[1][0],weights[2]*vectors[2][0]]
      ==[F(1,3),F(1,6),F(1,6)], ['1/3','1/6','1/6'])

# Rational coefficient polynomials in A,B: represent monomials by exponents.
def padd(x,y):
    keys=x.keys()|y.keys()
    return {k:x.get(k,F(0))+y.get(k,F(0)) for k in keys
            if x.get(k,F(0))+y.get(k,F(0))}
def pscale(x,c):
    return {k:v*c for k,v in x.items() if v*c}
def pmul(x,y):
    r={}
    for k,v in x.items():
        for l,w in y.items():
            m=tuple(a+b for a,b in zip(k,l))
            r[m]=r.get(m,F(0))+v*w
    return {k:v for k,v in r.items() if v}
A={(1,0):F(1)}; B={(0,1):F(1)}
target2=padd(pscale(pmul(A,A),F(3)),pmul(B,B))
case1=padd(pmul(pscale(A,F(2)),pscale(A,F(2))),pscale(target2,F(-1)))
case2=padd(pmul(padd(A,B),padd(A,B)),pscale(target2,F(-1)))
check('directional case A >= B >= 0 factors',
      case1==pmul(padd(A,pscale(B,F(-1))),padd(A,B)),
      '(2A)^2-(3A^2+B^2)=(A-B)(A+B) >= 0')
check('directional case B >= A >= 0 factors',
      case2==pscale(pmul(A,padd(B,pscale(A,F(-1)))),F(2)),
      '(A+B)^2-(3A^2+B^2)=2A(B-A) >= 0')

# Field elements (r,s) stand for r+s*sqrt(3).
def q3(r=0,s=0): return (F(r),F(s))
def qadd(x,y): return (x[0]+y[0],x[1]+y[1])
def qmul(x,y): return (x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
transform=[q3(1),q3(0,F(1,3))]
iso_vectors=[(q3(F(4,9)),q3()),
             (q3(F(4,3)),q3(0,F(4,3))),
             (q3(F(-4,3)),q3(0,F(4,3)))]
transformed=[]
for v in iso_vectors:
    w=[qmul(transform[i],v[i]) for i in range(2)]
    cv=[[qmul(w[i],w[j]) for j in range(2)] for i in range(2)]
    if not all(x[1]==0 for row in cv for x in row):
        raise ValueError('Coordinate transform retained irrational coefficients')
    transformed.append([[x[0] for x in row] for row in cv])
new_target=[[qmul(qmul(transform[i],q3(F(1,3) if i==j else 0)),transform[j])
             for j in range(2)] for i in range(2)]
check('coordinate transform diag(1,1/sqrt(3)) maps every covariance exactly',
      transformed==components and new_target==[[q3(v) for v in row] for row in sigma],
      'Diagonal component 2 vector changes sign without changing its outer product.')

t=F(1,24)
lower=t/F(12)-F(27,32)*t*t
delta=F(data['positive_definite_corollary']['delta'])
smooth_lower=lower-F(6,10000)
check('main rational lower bound', lower==F(data['expectation_gap_lower_bound'])>0, str(lower))
check('smoothed rational lower bound',
      smooth_lower==F(data['positive_definite_corollary']['expectation_gap_lower_bound'])>0,
      str(smooth_lower))
check('specified square root of smoothing variance', delta==F(1,10000)**2, str(delta))
for i,cov in enumerate([sigma]+components):
    sm=[[cov[j][k]+(delta if j==k else 0) for k in range(2)] for j in range(2)]
    check(f'smoothed covariance {i} positive definite by principal minors',
          sm[0][0]>0 and det(sm)>0, {'first_minor':str(sm[0][0]),'determinant':str(det(sm))})
affine_forms=[[F(x) for x in form] for form in data['witness_affine_slopes_and_intercepts']]
check('six affine slopes have Euclidean norms at most two',
      len(affine_forms)==6 and max(a*a+b*b for a,b,c in affine_forms)==4,
      {'squared_norms':[str(a*a+b*b) for a,b,c in affine_forms]})

# A further sanity check: covariance comparison is necessary, and holds here.
mixture_cov=[[sum(p*cov[i][j] for p,cov in zip(weights,components)) for j in range(2)] for i in range(2)]
covdiff=[[mixture_cov[i][j]-sigma[i][j] for j in range(2)] for i in range(2)]
check('mixture covariance dominates target', covdiff[0][0]>0 and det(covdiff)>0,
      {'mixture_covariance':[[str(x) for x in row] for row in mixture_cov],
       'difference':[[str(x) for x in row] for row in covdiff]})

# X=(R cos(theta)/sqrt(3),R sin(theta)/3), R Rayleigh, theta uniform.
# E[(k R-t)_+]=k*sqrt(pi/2)*erfc(t/(sqrt(2)*k)) for k>0.
# Symmetry reduces the target correction to theta in [0,pi/3].
c=math.sqrt(2/math.pi)
threshold=float(t)
def target_angular(threshold):
    def integrand(theta):
        k=math.cos(theta)/math.sqrt(3)-math.sin(theta)/3
        return 0.0 if k<=0 else k*erfc(threshold/(math.sqrt(2)*k))
    integral, err, info=quad(integrand,0,math.pi/3,epsabs=2e-14,epsrel=2e-14,full_output=1)
    return 2*c/3+c*integral,c*err,info['neval']

ex,ex_err,neval=target_angular(threshold)
# Independently integrate the mixture's horizontal excess directly over Z>=0.
horizontal,hor_err,hor_info=quad(lambda z:(4*z/9-threshold)*c*math.exp(-z*z/2),
                              threshold*9/4,math.inf,epsabs=2e-14,epsrel=2e-14,full_output=1)
ey=3*horizontal/4+2*c/3
ex0,ex0_err,ex0_n=target_angular(0)
report={'status':'independent_exact_checks_passed_optional_uncertified_numerics',
        'environment':{'python':platform.python_version(),'scipy':scipy.__version__},
        'exact_checks':checks,'exact_check_count':len(checks),
        'numerical_method':'Rayleigh radial identity and adaptive angular quadrature for X; direct half-normal integration for Y.',
        'target_expectation':ex,'mixture_expectation':ey,'gap':ex-ey,
        'estimated_target_quadrature_error':ex_err,'estimated_mixture_quadrature_error':3*hor_err/4,
        'integrand_evaluations':{'target':neval,'mixture_horizontal':hor_info['neval']},
        'zero_threshold_control':{'target':ex0,'known_common':c,'difference':ex0-c},
        'rational_lower_bound_decimal':float(lower),'smoothed_lower_bound_decimal':float(smooth_lower),
        'limitations':['Exact algebra checks do not constitute a formal proof of the probability arguments.',
                       'Numerics use double precision; quadrature error estimates are not certified.']}
print(json.dumps(report,indent=2))
