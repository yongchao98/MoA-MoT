import math

def solve_zeros_of_vector_field():
    """
    Calculates the least number of zeros a vector field can have on a compact manifold M.
    This is a theoretical problem. The user is expected to know the values of chi(M) and chi(dM).
    We will use some example values.
    """
    # Let's use the 3-ball (D^3) as an example.
    # M = D^3, the 3-dimensional ball.
    # The boundary of M is the 2-sphere, dM = S^2.
    chi_M = 1  # Euler characteristic of D^3
    chi_dM = 2 # Euler characteristic of S^2

    print(f"Let's consider an example: a 3-dimensional ball M = D^3.")
    print(f"The Euler characteristic of M is chi(M) = {chi_M}.")
    print(f"The boundary of M is the 2-sphere dM = S^2, with Euler characteristic chi(dM) = {chi_dM}.")

    # The least number of zeros a vector field can have on M is given by |chi(M)|.
    min_zeros = abs(chi_M)

    print("\nThe least number of zeros is given by the absolute value of the Euler characteristic of the manifold M, |chi(M)|.")
    print(f"Calculation: |{chi_M}| = {min_zeros}")
    print(f"So, the least number of zeros is {min_zeros}.")

solve_zeros_of_vector_field()