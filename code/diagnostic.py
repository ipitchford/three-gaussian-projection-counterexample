"""Independent diagnostic replay of the three-line asymmetric candidate.

Pure standard-library implementation; no numerical result is a proof.
No broader random search ran before this candidate was supplied.
Run: python3 code/diagnostic.py
Uses isotropic coordinates before the invertible change y -> y/sqrt(3).
"""

import json
import argparse
import math
import platform
import time
from pathlib import Path


ROOT_TWO_OVER_PI = math.sqrt(2 / math.pi)
THRESHOLD = 1 / 24


def half_normal_stoploss(sigma, threshold):
    """E[(|sigma Z|-threshold)_+] for threshold >= 0."""
    return (sigma * ROOT_TWO_OVER_PI * math.exp(-threshold**2 / (2 * sigma**2))
            - threshold * math.erfc(threshold / (math.sqrt(2) * sigma)))


def conditional_integrand(b, threshold=THRESHOLD):
    # X~N(0,I/3): A=|X_1| has sd 1/sqrt(3), B=|X_2|/sqrt(3) has sd 1/3.
    # max(A+B-t,2B)=2B+(A-B-t)_+.
    density_b = 3 * ROOT_TWO_OVER_PI * math.exp(-4.5 * b * b)
    return density_b * half_normal_stoploss(1 / math.sqrt(3), b + threshold)


def simpson(n, cutoff, threshold=THRESHOLD):
    if n <= 0 or n % 2:
        raise ValueError('Simpson interval count must be positive and even')
    step = cutoff / n
    odd = math.fsum(conditional_integrand(j * step, threshold) for j in range(1, n, 2))
    even = math.fsum(conditional_integrand(j * step, threshold) for j in range(2, n, 2))
    integral = step / 3 * (conditional_integrand(0, threshold)
                           + conditional_integrand(cutoff, threshold) + 4 * odd + 2 * even)
    return (2 / 3) * ROOT_TWO_OVER_PI + integral


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', action='store_true', help='Emit one JSON document only')
    args = parser.parse_args()
    started = time.monotonic()
    p = [3 / 4, 1 / 8, 1 / 8]
    sd = [4 / 9, 8 / 3, 8 / 3]
    directions = [[1.0, 0.0], [0.5, math.sqrt(3) / 2], [-0.5, math.sqrt(3) / 2]]
    covariances = [[[sd[i]**2 * directions[i][j] * directions[i][k]
                     for k in range(2)] for j in range(2)] for i in range(3)]
    ey = p[0] * half_normal_stoploss(sd[0], THRESHOLD) + (p[1] + p[2]) * sd[1] * ROOT_TWO_OVER_PI
    rows = []
    for cutoff in (2.0, 3.0, 4.0):
        for intervals in (512, 2048, 8192, 32768):
            ex = simpson(intervals, cutoff)
            rows.append({'cutoff_b': cutoff, 'intervals': intervals,
                         'target_expectation': ex, 'mixture_expectation': ey,
                         'gap_target_minus_mixture': ex - ey})
    # The integrand tail is bounded by E[A] P(B>cutoff); this is analytic,
    # but no rigorous floating-point evaluation of this bound is claimed.
    tail_bounds = {str(cutoff): ROOT_TWO_OVER_PI / math.sqrt(3) * math.erfc(3 * cutoff / math.sqrt(2))
                   for cutoff in (2.0, 3.0, 4.0)}
    direction_checks = []
    for i in range(3):
        total = math.fsum(p[j] * sd[j] * abs(directions[i][0] * directions[j][1]
                                           - directions[i][1] * directions[j][0]) for j in range(3))
        direction_checks.append({'perpendicular_to_component': i,
                                 'mixture_directional_sd': total,
                                 'target_directional_sd': 1 / math.sqrt(3),
                                 'margin': total - 1 / math.sqrt(3)})
    target_zero = simpson(32768, 4, 0)
    target_gap = rows[-1]['gap_target_minus_mixture']
    result = {
        'status': 'uncertified_numerical_diagnostic',
        'broader_search_status': 'No broad random search executed; redirected to supplied analytic candidate before implementation.',
        'target_covariance': [[1/3, 0], [0, 1/3]],
        'probabilities': p, 'component_sds': sd, 'directions': directions,
        'component_covariances': covariances,
        'witness': 'max(abs(x)+abs(y)/sqrt(3)-1/24, 2*abs(y)/sqrt(3))',
        'convexity': 'Maximum of the six affine forms +/-x +/-y/sqrt(3)-1/24 and +/-2y/sqrt(3).',
        'method': 'Conditional half-normal stop-loss formula and composite Simpson on B>=0; math.fsum summation.',
        'quadrature': rows, 'analytic_tail_bound_float_evaluations': tail_bounds,
        'rank_one_direction_tests': direction_checks,
        'zero_threshold_control': {'target_expectation': target_zero,
                                   'known_common_expectation': ROOT_TWO_OVER_PI,
                                   'difference': target_zero - ROOT_TWO_OVER_PI},
        'supplied_analytic_lower_bound': 37 / 18432,
        'numerical_gap_minus_supplied_lower_bound': target_gap - 37 / 18432,
        'limitations': [
            'No certified Simpson discretization or floating-point error bounds.',
            'Directional calculations use floating-point square roots; the three-check reduction needs its analytic proof.',
            'The supplied rational lower bound is compared numerically, not independently proved by this script.'
        ],
        'python': platform.python_version(), 'seconds': time.monotonic() - started
    }
    print(json.dumps(result, indent=2))
    if args.json:
        return
    print(json.dumps({'status': result['status'], 'target_expectation': rows[-1]['target_expectation'],
                      'mixture_expectation': ey, 'gap': target_gap,
                      'supplied_bound': 37 / 18432, 'seconds': result['seconds']}, indent=2))


if __name__ == '__main__':
    main()
