"""Finite internal checks only; the continuum/probability proof is in paper.tex."""
from pathlib import Path
from fractions import Fraction as F
import json, subprocess, sys, tempfile
R=Path(__file__).resolve().parents[1]
flags=['-O'] if sys.flags.optimize else []
def run(script,*args,ok=True):
    p=subprocess.run([sys.executable,*flags,str(R/'code'/script),*map(str,args)],capture_output=True,text=True,timeout=120)
    if (p.returncode==0)!=ok: raise RuntimeError(f'Unexpected exit for {script}: {p.returncode}\n{p.stderr}')
    return p.stdout
d=json.loads((R/'data/counterexample.json').read_text())
expected=[[str(a),str(b),'-1/24'] for a in (1,-1) for b in (1,-1)]+[['0','2','0'],['0','-2','0']]
if d['witness_affine_slopes_and_intercepts']!=expected: raise ValueError('Witness differs from manuscript')
if d['target_covariance']!=[['1/3','0'],['0','1/9']]: raise ValueError('Target differs from manuscript')
ind=json.loads(run('reviewer_checks.py','--submission',R))
num=json.loads(run('diagnostic.py','--json'))
if len(num['quadrature'])!=12 or abs(num['quadrature'][-1]['gap_target_minus_mixture']-ind['gap'])>1e-12:
    raise ValueError('Two numerical diagnostics disagree')
if abs(num['zero_threshold_control']['difference'])>1e-12: raise ValueError('Zero threshold control failed')
with tempfile.TemporaryDirectory(prefix='three-gaussian-negative-') as td:
    root=Path(td);(root/'data').mkdir()
    for name,mutate in [('wrong-weight',lambda x:x['weights'].__setitem__(0,'2/3')),('wrong-gap',lambda x:x.__setitem__('expectation_gap_lower_bound','1/10'))]:
        x=json.loads(json.dumps(d));mutate(x);(root/'data/counterexample.json').write_text(json.dumps(x))
        run('reviewer_checks.py','--submission',root,ok=False)
        print('Rejected semantic mutation:',name)
# Exact family sample checks supplement, but do not replace, the analytic formula.
for p in (F(7,10),F(3,4),F(9,10)):
    q=(1-p)/3;r=1-p-q
    if [p/(3*p),q/(6*q),r/(6*r)]!=[F(1,3),F(1,6),F(1,6)]: raise ValueError('Family coefficient error')
    t=(p-F(2,3))/(3*p*p)
    if (p-F(2,3))*t-F(3,2)*p*p*t*t<=0: raise ValueError('Family strictness error')
print(json.dumps({'status':'passed','exactChecks':ind['exact_check_count'],'numericalRows':12,'negativeControls':2,'optimized':bool(sys.flags.optimize),'gapDiagnostic':ind['gap'],'scope':'Internal finite checks; no formal proof or certified quadrature.'}))
