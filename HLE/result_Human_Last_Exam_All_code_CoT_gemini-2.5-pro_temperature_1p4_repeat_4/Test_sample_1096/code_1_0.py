import sympy as sp

def solve_optimal_step_sizes():
    """
    This function derives and prints the optimal step sizes (gamma_1, gamma_2)
    for a two-step gradient descent method on k-smooth and 1-strongly convex functions.
    The goal is to minimize the worst-case convergence ratio ||x_2 - x*||^2 / ||x_0 - x*||^2.
    """
    # M is the smoothness constant, m is the strong convexity constant.
    # kappa is the condition number M/m.
    # The problem assumes m=1 and M=kappa.
    kappa = sp.Symbol('kappa', positive=True, real=True)

    # The optimal step sizes are derived from the roots of the optimal polynomial,
    # which is a scaled Chebyshev polynomial of degree 2.
    # The derivation leads to the following expressions for gamma_1 and gamma_2.
    
    # Numerator of the expression for the step sizes
    # 4*(kappa + 1) - 2*sqrt(2)*(kappa - 1)
    # 4*(kappa + 1) + 2*sqrt(2)*(kappa - 1)
    
    # We can write out each part of the formula to satisfy the "output each number" requirement
    four = sp.Integer(4)
    one = sp.Integer(1)
    two = sp.Integer(2)
    six = sp.Integer(6)
    
    numerator_part1 = four * (kappa + one)
    numerator_part2 = two * sp.sqrt(two) * (kappa - one)
    
    denominator = kappa**two + six * kappa + one

    # The two optimal step sizes
    gamma_1 = (numerator_part1 - numerator_part2) / denominator
    gamma_2 = (numerator_part1 + numerator_part2) / denominator

    print("The best choice for the pair (gamma_1, gamma_2) is given by the following formulas:")
    print("-" * 20)
    
    print("gamma_1 = ")
    sp.pprint(gamma_1)
    
    print("\ngamma_2 = ")
    sp.pprint(gamma_2)
    
    print("-" * 20)
    print("Note: The order of gamma_1 and gamma_2 can be swapped.")
    print("kappa is the condition number M/m.")

solve_optimal_step_sizes()