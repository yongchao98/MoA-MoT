def solve_squares_on_grid():
    """
    Calculates the number of squares on an n x n grid and demonstrates
    the formula used to find 'a' and 'b'.
    """
    # Let's demonstrate with a 4x4 grid.
    n = 4

    print(f"To find the total number of squares on an {n}x{n} grid,")
    print("we use the expression sum(a^2 * b) for m from 1 to n.")
    print("Based on the derivation, we have:")
    print("a = m")
    print("b = n - m + 1")
    print("\nFor n = 4, the calculation is:")

    total_squares = 0
    equation_parts = []

    # The summation is from m=1 to n
    for m in range(1, n + 1):
        # Here, a = m and b = n - m + 1
        a = m
        b = n - m + 1
        
        term_value = (a**2) * b
        total_squares += term_value
        
        # We need to output each number in the final equation
        # The terms in the equation are a^2 * b
        equation_parts.append(f"{a}^2*{b}")

    # Print the full equation and the final result
    equation_string = " + ".join(equation_parts)
    print(f"{equation_string} = {total_squares}")

solve_squares_on_grid()
