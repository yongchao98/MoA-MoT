import math

def solve_zeros_on_manifold():
    """
    Calculates the least number of zeros a vector field can have on a
    compact manifold M with a non-empty boundary.
    """
    # We use the 3-dimensional disk (M = D^3) as a concrete example.
    # The Euler characteristic of the 3-disk is 1.
    chi_M = 1

    # The boundary of D^3 is the 2-sphere (dM = S^2).
    # The Euler characteristic of the 2-sphere is 2.
    # This value is included for context but is not used in the final calculation.
    chi_dM = 2

    # The least number of zeros a vector field on M can have is the
    # absolute value of the Euler characteristic of M, |chi(M)|.
    # This result comes from the Poincaré-Hopf theorem for manifolds with boundary.
    min_zeros = abs(chi_M)

    print("The problem is to find the least number of zeros a vector field can have on a compact manifold M with a non-empty boundary.")
    print("The answer is given by the absolute value of the Euler characteristic of M, |chi(M)|.")
    print("\n--- Example Calculation ---")
    print(f"Let M be the 3-disk (D^3).")
    print(f"The Euler characteristic of M is chi(M) = {chi_M}.")
    print(f"The Euler characteristic of its boundary (S^2) is chi(dM) = {chi_dM}.")
    print("\nThe final equation for the minimum number of zeros is:")
    print("min_zeros = |chi(M)|")
    print("\nPlugging in the numbers for our example:")
    # The user request was to output each number in the final equation.
    print(f"min_zeros = |{chi_M}| = {min_zeros}")

solve_zeros_on_manifold()