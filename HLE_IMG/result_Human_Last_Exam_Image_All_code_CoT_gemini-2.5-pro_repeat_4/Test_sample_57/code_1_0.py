def solve():
    """
    This function provides the formula for the number of ways to fill the shape
    with 2x1 rectangles.
    """
    # The number of ways to tile the given L-shape of size n is given by the formula:
    # a_n = 2 * F_{2n-2}, where F_k is the k-th Fibonacci number.
    
    # The constant coefficients and values in the formula are:
    coefficient = 2
    n_coefficient = 2
    constant_term = -2
    
    # Print the final formula, showing each numerical component.
    print("The number of ways to fill out the shape is expressed by the formula:")
    print(f"{coefficient} * F_{{{n_coefficient}*n + ({constant_term})}}")
    print("\nWhich simplifies to:")
    print(f"{coefficient} * F_{{{n_coefficient}*n - 2}}")

solve()