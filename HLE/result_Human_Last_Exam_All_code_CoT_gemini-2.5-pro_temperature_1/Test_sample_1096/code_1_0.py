import sympy

def find_optimal_steps():
    """
    This function computes and prints the symbolic formulas for the optimal
    step sizes (gamma_1, gamma_2) for a two-step gradient descent method.

    The problem is to minimize the worst-case convergence rate for M-smooth
    and mu-strongly convex functions. We use the settings mu=1 and M=kappa.
    """

    # Define kappa as a symbol for the condition number M/mu.
    # The problem specifies mu=1, so M=kappa.
    kappa = sympy.Symbol('kappa')
    mu = 1

    # The optimal step sizes gamma_1 and gamma_2 are the reciprocals of the
    # roots of the scaled Chebyshev polynomial of degree 2.
    # The roots r_1, r_2 of the polynomial P(lambda) are given by:
    # r_{1,2} = (M+mu)/2 \pm (M-mu)/(2*sqrt(2))
    # We substitute M=kappa and mu=1.
    # The step sizes are gamma_1 = 1/r_1 and gamma_2 = 1/r_2.

    sqrt2 = sympy.sqrt(2)

    # Simplified formulas for gamma_1 and gamma_2 derived from the roots.
    # These are two distinct values that form the optimal pair.
    # The order does not matter.
    gamma1_expr = (2 * sqrt2) / ((sqrt2 + 1) * kappa + (sqrt2 - 1) * mu)
    gamma2_expr = (2 * sqrt2) / ((sqrt2 - 1) * kappa + (sqrt2 + 1) * mu)

    print("The optimal choice for the pair of step sizes (gamma_1, gamma_2) is given by the following two formulas:")
    
    print("\nFormula 1:")
    sympy.pprint(gamma1_expr, use_unicode=True)
    
    print("\nFormula 2:")
    sympy.pprint(gamma2_expr, use_unicode=True)

    print("\nNote: The variable 'kappa' represents the condition number M (since mu=1).")
    print("Each number in the final equations is printed above.")

if __name__ == '__main__':
    find_optimal_steps()