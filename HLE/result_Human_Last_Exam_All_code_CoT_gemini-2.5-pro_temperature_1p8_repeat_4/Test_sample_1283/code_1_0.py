def solve_max_solutions():
    """
    Calculates the maximum number of solutions for the given problem.
    """
    
    # (a) The general formula is derived as d_P + d_Q + 2.
    # We will represent this as a string.
    general_formula = "d_P + d_Q + 2"
    
    # (b) We are given specific degrees for the polynomials P(x) and Q(x).
    d_P = 3
    d_Q = 2
    
    # Calculate the maximum number of solutions for the specific case.
    max_solutions = d_P + d_Q + 2
    
    # Print the results in the required format.
    # First, the expression for part (a).
    print(f"(a) {general_formula}")
    
    # Second, the detailed calculation and result for part (b).
    # The instructions require printing each number in the final equation.
    print(f"(b) {d_P} + {d_Q} + 2 = {max_solutions}")

solve_max_solutions()