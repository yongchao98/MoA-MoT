import math

def solve_optimal_steps():
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for a given condition number kappa.
    """
    try:
        # Ask user for the condition number kappa (M), assuming mu=1.
        kappa_str = input("Please enter the condition number kappa (where kappa = M and mu = 1): ")
        kappa = float(kappa_str)
        if kappa <= 1:
            print("Error: The condition number kappa must be greater than 1.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a numerical value for kappa.")
        return

    # M is the smoothness constant, mu is the strong convexity constant.
    M = kappa
    mu = 1.0

    # The problem is to find gamma_1 and gamma_2 to minimize the convergence rate.
    # The analysis leads to a minimax problem solved by Chebyshev polynomials.
    # The optimal step sizes are the reciprocals of the roots of the scaled Chebyshev polynomial.
    # The formulas for the optimal steps are:
    # gamma_1 = 2 / ((M + mu) + (M - mu) / sqrt(2))
    # gamma_2 = 2 / ((M + mu) - (M - mu) / sqrt(2))
    # Note: The order of gamma_1 and gamma_2 can be swapped.

    sqrt2 = math.sqrt(2)
    
    denominator_1 = (M + mu) + (M - mu) / sqrt2
    denominator_2 = (M + mu) - (M - mu) / sqrt2

    gamma_1 = 2.0 / denominator_1
    gamma_2 = 2.0 / denominator_2

    print(f"\nFor a function with M = {M} and mu = {mu}:")
    print("The best choice for the pair (gamma_1, gamma_2) is:")
    # The final equation is the value of gamma1 and gamma2. We print them number by number.
    print(f"gamma_1 = {gamma_1}")
    print(f"gamma_2 = {gamma_2}")
    
solve_optimal_steps()