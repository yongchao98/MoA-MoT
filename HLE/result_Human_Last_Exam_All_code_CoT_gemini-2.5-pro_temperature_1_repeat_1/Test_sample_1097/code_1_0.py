import numpy as np
from scipy.optimize import fsolve

def solve_and_explain():
    """
    This function derives the set M for the Banach Fixed-Point Theorem application,
    calculates the specific numerical constraints, and prints the full explanation.
    """

    # Derive the numerical constraints for R
    # The two conditions for R are:
    # 1. Contraction: exp(R) / 8 < 1  => R < ln(8)
    # 2. Invariance (T(M) subset M): exp(R) / 8 <= R

    # We find the lower bound for R by solving the equality R = exp(R) / 8.
    # We define the function g(R) = R - exp(R) / 8 and find its root.
    def g(R):
        return R - np.exp(R) / 8

    # An initial guess of 0.1 is reasonable as g(0) is negative and g(1) is positive.
    R_lower_bound = fsolve(g, 0.1)[0]

    # The upper bound is ln(8).
    R_upper_bound = np.log(8)

    # Print the detailed explanation and the final answer.
    print("""
To prove the existence and uniqueness of a global solution to the boundary value problem
    u''(x) - exp(u(x)) = 0, x ∈ (0, 1)
    u(0) = u(1) = 0
using the Banach Fixed-Point Theorem, you should define the set M as follows:

1.  **The Banach Space (X, ||·||_∞):**
    First, we define the space X as the set of all continuous functions on the interval [0, 1] that are zero at the endpoints x=0 and x=1.
    X = C₀([0, 1]) = {u ∈ C([0, 1]) | u(0) = u(1) = 0}
    This is a complete metric space (a Banach space) when equipped with the supremum norm:
    ||u||_∞ = max_{x ∈ [0, 1]} |u(x)|

2.  **The Fixed-Point Operator T:**
    The boundary value problem is equivalent to the fixed-point integral equation u = Tu, where T is the operator:
    (Tu)(x) = ∫[0,1] G(x, s) * exp(u(s)) ds
    Here, G(x, s) is the Green's function for u'' with zero boundary conditions, which is G(x,s) = s(x-1) for s ≤ x and G(x,s) = x(s-1) for s ≥ x.

3.  **The Set M:**
    To ensure the conditions of the Banach theorem are met, we define M as a closed ball of radius R > 0 within the space X. The set M is defined as:

    M = {u ∈ C₀([0, 1]) | ||u||_∞ ≤ R}

    where the radius R must be chosen such that the operator T maps M into itself (T(M) ⊆ M) and is a contraction on M. These requirements lead to the following conditions on R:
    
    a) For T to be a contraction: (1/8) * exp(R) < 1  =>  R < ln(8)
    b) For T to map M to itself: (1/8) * exp(R) ≤ R

    Combining these, we need to find a value R that satisfies R ≥ exp(R)/8 and R < ln(8).
    The valid range for R is [R₀, ln(8)), where R₀ is the smallest positive root of the equation R = exp(R)/8.

    Numerically, the bounds for R are:
    Lower bound, R₀ ≈ {r_lower:.5f}
    Upper bound, ln(8) ≈ {r_upper:.5f}

    Therefore, a valid set M is any closed ball in C₀([0, 1]) with a radius R chosen from the interval [{r_lower:.5f}, {r_upper:.5f}).

    For example, choosing R=1 is a valid choice since 1 is in this interval.
    """.format(r_lower=R_lower_bound, r_upper=R_upper_bound))

if __name__ == '__main__':
    solve_and_explain()
    print("<<<M = {u ∈ C₀([0, 1]) | ||u||_∞ ≤ R}, where R is any constant in the interval [{:.5f}, {:.5f})>>>".format(fsolve(lambda R: R - np.exp(R)/8, 0.1)[0], np.log(8)))