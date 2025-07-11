import sympy

def solve_limit_problem():
    """
    This function calculates the limit of n*P(n) as n -> infinity based on the
    analytical derivation explained above. It uses the sympy library for
    symbolic mathematics to compute the final limit.
    """

    # Define k as a symbolic variable, representing a positive value tending to infinity.
    k = sympy.Symbol('k', positive=True)

    # n is defined as 6k. As n -> oo, k -> oo.
    n = 6 * k

    # The expression for P(n) for large n (and k) is derived from the
    # Central Limit Theorem.
    # P(n) is approximately 1 - exp(-1/(3*k)).
    P_n_approx = 1 - sympy.exp(-1 / (3 * k))

    # We are interested in the limit of the product n * P(n).
    expression_to_limit = n * P_n_approx

    print("Based on the Central Limit Theorem, the expression for n*P(n) for large n=6k is:")
    
    # The following print statement outputs the symbolic expression for n*P(n).
    # The equation for the limit is: L = lim_{k->inf} [expression]
    # This prints the expression part, which contains all the relevant numbers (6, 1, 3)
    # from the derived formula.
    sympy.pprint(expression_to_limit, use_unicode=True)

    # Calculate the limit of the expression as k goes to infinity.
    limit_value = sympy.limit(expression_to_limit, k, sympy.oo)

    print("\nThe limit of the expression as k -> oo is:")
    print(limit_value)

solve_limit_problem()