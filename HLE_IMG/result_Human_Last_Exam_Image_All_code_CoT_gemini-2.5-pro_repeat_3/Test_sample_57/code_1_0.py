def solve_tiling():
    """
    This function provides the formula for the number of ways to tile the given shape.

    The shape is an n x n square with an (n-2) x (n-2) square removed from a corner.
    The number of ways to tile this shape with 2x1 dominoes is a known result
    in combinatorics, given by the formula 2 * (F_n)^2, where F_n is the n-th
    Fibonacci number.
    """

    # The formula involves the following numbers:
    coefficient = 2
    exponent = 2

    # We print the final formula.
    print("The number of ways to fill the shape is given by the formula:")
    print(f"{coefficient} * (F_n)^{exponent}")

solve_tiling()
