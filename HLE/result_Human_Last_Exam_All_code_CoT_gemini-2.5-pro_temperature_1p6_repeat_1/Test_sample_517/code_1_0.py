import sympy

def solve_probability_limit():
    """
    This function demonstrates the final step of the limit calculation symbolically.

    Let a(x) be the potential kernel of the 2D simple random walk.
    Let x0 = (0, 1) be the starting point.
    Let Tn = (n, 0) be the target point for large n.
    
    The probability p_n is approximated by:
    p_n ≈ (a(x0) - a(x0 - Tn) + a(Tn)) / (2 * a(x0))
    
    We need to find the limit of p_n as n -> ∞.
    """

    # Define symbols for the terms in the equation.
    # a_x0 represents a((0,1)), which is a fixed constant.
    a_x0 = sympy.Symbol('a_x0')
    
    # a_Tn represents a((n,0)). As n -> ∞, this term -> ∞.
    a_Tn = sympy.Symbol('a_Tn')
    
    # a_x0_minus_Tn represents a((0,1) - (n,0)) = a((-n,1)).
    # As n -> ∞, ||(-n,1)|| ≈ n, so a((-n,1)) ≈ a((n,0)).
    # In the limit, a((-n,1)) and a((n,0)) become equal.
    a_x0_minus_Tn = sympy.Symbol('a_x0_minus_Tn')

    # Expression for the numerator of p_n * a(x0)
    numerator = (a_x0 - a_x0_minus_Tn + a_Tn) / 2
    
    # Full expression for p_n
    p_n_expr = numerator / a_x0

    print("The expression for p_n is: " + str(p_n_expr))
    
    # In the limit n -> ∞, a(x0 - Tn) tends to the same value as a(Tn).
    # We can substitute a_x0_minus_Tn with a_Tn to find the limit.
    limit_expr = p_n_expr.subs(a_x0_minus_Tn, a_Tn)
    
    print("\nIn the limit as n -> infinity, a(x0 - Tn) and a(Tn) are equal.")
    print("Substituting a(x0 - Tn) with a(Tn), the expression becomes: " + str(limit_expr))
    
    # The final equation is derived from the simplification of the limit expression.
    final_eq_lhs = limit_expr.simplify()
    
    # The terms cancel out nicely.
    # We can represent the final simplified equation with its numbers.
    final_numerator = 1
    final_denominator = 2
    
    print(f"\nThe simplified limit is the fraction: {final_numerator} / {final_denominator}")

    final_answer = final_numerator / final_denominator
    print("\nFinal Result:")
    print(f"The final simplified equation is {final_numerator}/{final_denominator}")
    
solve_probability_limit()
