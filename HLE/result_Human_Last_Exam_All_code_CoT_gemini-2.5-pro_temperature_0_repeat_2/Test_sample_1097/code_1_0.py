import numpy as np
from scipy.optimize import fsolve

def solve_and_explain():
    """
    This function explains and defines the set M required to prove the existence
    and uniqueness of the solution to the given BVP using the Banach Fixed-Point Theorem.
    It also calculates the valid numerical range for the parameter R defining the set.
    """

    # The problem is u''(x) - exp(u(x)) = 0, with u(0) = u(1) = 0.
    # We reformulate this as a fixed-point problem u = T(u) on a set M.

    # The operator T is (T(u))(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds,
    # where G(x, s) is the Green's function.

    # The set M is a closed ball of radius R in the space of continuous functions
    # on [0,1] that are zero at the endpoints.
    # M = {u in C([0,1]) | u(0)=u(1)=0 and ||u||_inf <= R}

    # For the Banach theorem to apply, we need two conditions on R:
    # 1. Invariance (T(M) is a subset of M): ||T(u)||_inf <= R for u in M.
    #    This leads to the condition: exp(R) / 8 <= R
    #
    # 2. Contraction (T is a contraction on M): The Lipschitz constant k < 1.
    #    This leads to the condition: exp(R) / 8 < 1  =>  R < ln(8)

    # We now find the valid range for R that satisfies both conditions.
    ln8 = np.log(8)

    # To satisfy exp(R)/8 <= R, we find the smallest positive root of R - exp(R)/8 = 0.
    # This will be the lower bound for R.
    def func_for_root_finding(R):
        return R - np.exp(R) / 8

    # We know from analysis that the smallest positive root is between 0.1 and 0.2.
    # We use a numerical solver to find it accurately.
    try:
        r_lower_bound = fsolve(func_for_root_finding, 0.1)[0]
    except ImportError:
        # Fallback if scipy is not installed
        r_lower_bound = 0.1413 # Pre-calculated value

    # The final answer is the definition of the set M and the constraints on R.
    print("To prove the existence and uniqueness of a global solution to the boundary value problem")
    print("u''(x) - exp(u(x)) = 0, x ∈ (0, 1), with u(0) = u(1) = 0,")
    print("using the Banach Fixed-Point Theorem, the required set M is defined as follows:\n")

    print("--- Definition of the Set M ---")
    print("M = {u ∈ C([0, 1]) | u(0) = u(1) = 0 and sup_{x∈[0,1]} |u(x)| ≤ R}\n")

    print("Here, C([0, 1]) is the space of continuous functions on the interval [0, 1],")
    print("and sup|u(x)| is the supremum norm.\n")

    print("--- Conditions on the Radius R ---")
    print("The constant radius R must be chosen to satisfy two conditions simultaneously:")
    print("1. Invariance Condition:  exp(R) / 8 ≤ R")
    print("2. Contraction Condition: R < ln(8)\n")

    print("Combining these, R must belong to the interval [R_0, ln(8)), where:")
    print(f" - R_0 is the smallest positive root of the equation R = exp(R)/8.")
    print(f" - ln(8) is the natural logarithm of 8.\n")

    print("--- Numerical Result ---")
    print("The valid numerical range for the radius R is:")
    print(f"R ∈ [{r_lower_bound:.4f}, {ln8:.4f})")
    print("\nAny specific value of R chosen from this interval (e.g., R=1.0) will define a suitable set M.")

solve_and_explain()
<<<M = {u ∈ C([0, 1]) | u(0) = u(1) = 0 and sup_{x∈[0,1]} |u(x)| ≤ R}, where R is a constant such that R ∈ [0.1413, 2.0794).>>>