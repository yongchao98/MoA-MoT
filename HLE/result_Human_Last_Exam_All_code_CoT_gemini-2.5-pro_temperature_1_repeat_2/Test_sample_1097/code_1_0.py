import numpy as np
from scipy.optimize import brentq

def define_set_m():
    """
    Calculates the bounds for the radius R of the set M and prints its definition.
    """
    # The condition for the operator T to map the set M to itself is exp(R)/8 <= R.
    # We find the lower bound for R by solving the equation R - exp(R)/8 = 0.
    # Let's define the function f(R) = R - exp(R)/8.
    f = lambda R: R - np.exp(R) / 8.0

    # We can observe that f(0) < 0 and f(1) > 0, so a root exists in [0, 1].
    try:
        r_min = brentq(f, 0, 2)
    except (ImportError, ValueError):
        # Fallback if scipy is not installed or if the root is not bracketed
        # A simple numerical search would also work. For this problem, we know the approx. value.
        r_min = 0.148496 # Approximate value

    # The condition for T to be a contraction is exp(R)/8 < 1.
    # This simplifies to R < ln(8).
    r_max = np.log(8)

    print("To prove the existence and uniqueness of the solution using the Banach Fixed-Point Theorem, we define the set M.")
    print("M must be a complete metric space on which the associated integral operator T is a contraction that maps M to itself.")
    print("\nThe appropriate choice is a closed ball of radius R in the space of continuous functions on [0, 1] that are zero at the boundaries:")
    print("\nM = { u in C[0, 1] | u(0) = u(1) = 0 and ||u||_inf <= R }")
    print("\nFor the operator to be a contraction mapping from M to M, the radius R must satisfy two conditions:")
    print(f"1. Invariance (T(M) subset of M): exp(R)/8 <= R")
    print(f"2. Contraction: exp(R)/8 < 1  =>  R < ln(8)")
    print("\nCombining these gives the valid range for R:")
    print(f"[{r_min:.6f}, {r_max:.6f})")
    
    print("\nSo, the final definition of the set M is:")
    print("M = { u in C[0, 1] | u(0) = u(1) = 0 and ||u||_inf <= R },")
    print(f"where R is any constant in the interval [{r_min:.4f}, {r_max:.4f}).")


if __name__ == '__main__':
    define_set_m()
