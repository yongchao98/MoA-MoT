import numpy as np

def solve_walk_problem():
    """
    This function demonstrates that for k=2, it's possible to find measures
    that allow for a recurrent controlled random walk in d>=3 dimensions.
    It does this by constructing two covariance matrices whose convex
    combination is the identity matrix.
    """
    d = 3  # Dimension, must be >= 3

    print(f"Demonstrating the construction for d = {d} and k = 2.\n")

    # Define the two covariance matrices
    sigma1_diag = np.ones(d)
    sigma1_diag[0] = d
    Sigma1 = np.diag(sigma1_diag)

    sigma2_diag = np.ones(d)
    sigma2_diag[0] = 1/d
    Sigma2 = np.diag(sigma2_diag)

    # Define the coefficients for the convex combination
    alpha1 = 1 / (d + 1)
    alpha2 = d / (d + 1)

    # Calculate the convex combination
    M = alpha1 * Sigma1 + alpha2 * Sigma2

    # Print the equation
    np.set_printoptions(precision=3)
    print("The equation is: alpha1 * Sigma1 + alpha2 * Sigma2 = M\n")
    print(f"With alpha1 = {alpha1:.3f} and alpha2 = {alpha2:.3f}\n")
    
    print("Sigma1:\n", Sigma1, "\n")
    print("Sigma2:\n", Sigma2, "\n")
    
    print(f"{alpha1:.3f} * \n{Sigma1}\n+ {alpha2:.3f} * \n{Sigma2}\n= \n{M}")
    
    print("\nThe resulting matrix M is the identity matrix.")
    print("This shows that for k=2, we can find measures that allow for a recurrent walk.")
    print("Therefore, the maximal k for which the walk is always transient is 1.")

solve_walk_problem()