import sympy

def solve_limit_problem():
    """
    This function calculates the limit of n*P(n) as n -> infinity.
    
    The derivation shows that for large n, P(n) can be approximated by
    the expression 1 - exp(-2/n). This function uses this expression
    to compute the final limit symbolically.
    """
    
    # Define the symbol n for our limit calculation.
    n = sympy.Symbol('n', positive=True)
    
    # This is the expression for P(n) derived from the Central Limit Theorem.
    # The numbers in this expression are 1 and -2.
    P_n = 1 - sympy.exp(-2 / n)
    
    # This is the expression for which we need to find the limit.
    expression_to_limit = n * P_n
    
    # Calculate the limit as n approaches infinity.
    limit_value = sympy.limit(expression_to_limit, n, sympy.oo)
    
    # Print the final equation and its result.
    # The final equation is lim_{n->oo} n * (1 - exp(-2/n)) = 2.
    # The numbers involved are 1, -2, and the result 2.
    print("The expression for P(n) for large n is:")
    print("P(n) =", P_n)
    print("\nWe are calculating the limit of the following expression as n -> oo:")
    print("Expression =", expression_to_limit)
    print("\nThe final equation and its result are:")
    # Using sympy.pretty_print for a more mathematical output format
    final_equation = sympy.Eq(sympy.Limit(expression_to_limit, n, sympy.oo), limit_value)
    sympy.pprint(final_equation)

solve_limit_problem()