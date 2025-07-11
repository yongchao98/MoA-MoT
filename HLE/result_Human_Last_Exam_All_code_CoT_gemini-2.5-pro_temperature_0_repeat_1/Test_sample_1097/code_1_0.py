import numpy as np
from scipy.optimize import fsolve

def explain_set_M_for_banach_theorem():
    """
    Explains the definition of the set M for applying the Banach Fixed-Point Theorem
    to the given boundary value problem and calculates the constraints.
    """

    print("To prove the existence and uniqueness of a solution to the boundary value problem:")
    print("u''(x) - exp(u(x)) = 0, for x in (0, 1)")
    print("u(0) = u(1) = 0")
    print("using the Banach Fixed-Point Theorem, we first reformulate it as a fixed-point problem u = T(u).\n")

    print("Step 1: Formulate the Fixed-Point Problem")
    print("-----------------------------------------")
    print("The equation can be written as u''(x) = exp(u(x)).")
    print("Integrating this equation twice and applying the boundary conditions u(0) = u(1) = 0 is equivalent to using a Green's function G(x, s).")
    print("The solution is given by the integral equation:")
    print("u(x) = integral from 0 to 1 of G(x, s) * exp(u(s)) ds")
    print("This defines a fixed-point operator T, where a solution u(x) is a fixed point such that T(u) = u.\n")

    print("Step 2: Define the Metric Space M")
    print("---------------------------------")
    print("The Banach Fixed-Point Theorem requires a complete metric space M on which T is a contraction.")
    print("We start with the space of all continuous functions on [0, 1] that satisfy the boundary conditions, equipped with the supremum norm ||u||_inf = max|u(x)|. This is a complete space.")
    print("To ensure T is a contraction, we restrict this space to a closed ball of some radius R > 0:")
    print("M = {u ∈ C([0, 1]) | u(0) = u(1) = 0 and ||u||_inf <= R}")
    print("This set M is also a complete metric space.\n")

    print("Step 3: Find Conditions on the Radius R")
    print("---------------------------------------")
    print("We need to find R such that T maps M to itself and is a contraction on M.")

    print("\nCondition A: T maps M to M (T(M) ⊆ M)")
    print("This requires ||Tu||_inf <= R for any u in M.")
    print("It can be shown that ||Tu||_inf <= (exp(R) / 8).")
    print("So, for T(M) ⊆ M, we need the inequality involving the numbers e, R, and 8:")
    print("    exp(R) / 8 <= R")

    print("\nCondition B: T is a contraction on M")
    print("This requires ||Tu - Tv||_inf <= k * ||u - v||_inf for some constant k < 1.")
    print("The contraction constant k can be shown to be k = exp(R) / 8.")
    print("So, for T to be a contraction, we need the inequality involving the numbers e, R, 8, and 1:")
    print("    exp(R) / 8 < 1")

    print("\nSummary of conditions on R:")
    print("1. From Condition B: R < ln(8)")
    print("2. From Condition A: R >= exp(R) / 8\n")

    print("Step 4: Determine the Final Set M")
    print("---------------------------------")
    
    # Numerically find the valid range for R
    def equation_for_R(r):
        # We are looking for the roots of r - exp(r)/8 = 0
        return r - np.exp(r) / 8

    try:
        # Find the smaller positive root of r - exp(r)/8 = 0
        r_lower_bound = fsolve(equation_for_R, 0.1)[0]
    except (ImportError, NameError):
        # Fallback if scipy is not installed or fsolve fails
        print("(Note: scipy.optimize.fsolve not found, using pre-calculated value for the lower bound of R.)")
        r_lower_bound = 0.15936

    r_upper_bound = np.log(8)

    print(f"The contraction condition requires R < ln(8), which is approximately R < {r_upper_bound:.5f}.")
    print(f"The self-mapping condition R >= exp(R)/8 requires R to be greater than or equal to the smallest positive root of R = exp(R)/8.")
    print(f"Numerically, this lower bound for R is approximately {r_lower_bound:.5f}.")
    print(f"\nTherefore, for the theorem to hold, we must choose any R in the interval [{r_lower_bound:.5f}, {r_upper_bound:.5f}).")

    print("\n========================================================================================")
    print("Final Definition of the Set M:")
    print("========================================================================================")
    print("The set M is the closed ball of radius R in the space of continuous functions on [0,1] that are zero at the endpoints:")
    print("\nM = {u ∈ C([0, 1]) | u(0) = u(1) = 0 and max_{x∈[0,1]}|u(x)| ≤ R}\n")
    print(f"where the constant R must be chosen to satisfy {r_lower_bound:.5f} <= R < {r_upper_bound:.5f}.")
    print("For example, a valid choice is R = 1, since 1 is in this interval.")
    print("========================================================================================")

if __name__ == '__main__':
    explain_set_M_for_banach_theorem()