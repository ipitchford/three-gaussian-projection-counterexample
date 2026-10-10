# Projection tests fail for three centred Gaussian components

Anonymous · 1.1.0-candidate · DOI 10.5281/zenodo.23284191

A written-proof candidate separates every projected convex comparison from full joint convex order for one Gaussian versus three centred Gaussian components. This is not a claim that all such mixtures fail, nor a classification of all weights. The two-versus-three boundary relies on Guéant's cited two-component theorem.

Start with [paper.pdf](paper.pdf), [AI_INDEX.md](AI_INDEX.md), [ASSURANCE.md](ASSURANCE.md) and [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md).

## Reproduce
Python 3.13: `python -m pip install -r requirements.txt`, then `python code/run_checks.py` and `python -O code/run_checks.py`. `python code/check_manifest.py` verifies shipped file hashes. `python code/diagnostic.py --json` emits a single JSON object. Build the paper with two passes of `pdflatex paper.tex`.

The numerical integrations are diagnostics, not certified bounds. The continuum-of-directions argument and strict rational lower bound are analytic in the paper. The original input and supplied review are preserved privately; historical process logs and third-party paper bodies are not redistributed as original research.
