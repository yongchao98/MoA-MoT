import math

def solve_optimal_gammas():
    """
    This function provides the formulas for the optimal pair of step sizes
    (gamma_1, gamma_2) for a two-step gradient descent algorithm.

    The optimization minimizes the worst-case convergence rate for M-smooth and
    mu-strongly convex functions, where M=kappa and mu=1.

    The solution is derived using Chebyshev polynomials of degree 2. The inverse
    of the step sizes (gamma_1, gamma_2) are the roots of the optimal polynomial.
    The roots are found to be:
    r_1 = (kappa+1)/2 + (kappa-1)/(2*sqrt(2))
    r_2 = (kappa+1)/2 - (kappa-1)/(2*sqrt(2))

    The optimal step sizes are the reciprocals of these roots.
    """

    # We express the formulas as strings for clarity.
    # The order of gamma_1 and gamma_2 can be swapped.
    # gamma_1 corresponds to the reciprocal of r_1
    gamma1_numerator = "2 * sqrt(2)"
    gamma1_denominator = "(sqrt(2) + 1) * kappa + (sqrt(2) - 1)"
    
    # gamma_2 corresponds to the reciprocal of r_2
    gamma2_numerator = "2 * sqrt(2)"
    gamma2_denominator = "(sqrt(2) - 1) * kappa + (sqrt(2) + 1)"

    print("The best choice for the pair of step sizes (gamma_1, gamma_2) is given by:")
    print("-" * 70)
    print(f"         {gamma1_numerator}")
    print(f"gamma_1 = {'-' * (len(gamma1_denominator) + 4)}")
    print(f"          {gamma1_denominator}")
    print("\n")
    print(f"         {gamma2_numerator}")
    print(f"gamma_2 = {'-' * (len(gamma2_denominator) + 4)}")
    print(f"          {gamma2_denominator}")
    print("-" * 70)
    print("\nNote on the notation 'S':")
    print("The term S = sqrt(kappa^2 + (kappa-1)^2) mentioned in the problem description")
    print("does not appear in the standard derivation for this problem and is likely unrelated.")
    print("The solution provided above is the established result from optimization theory.")

if __name__ == '__main__':
    solve_optimal_gammas()
