def solve_squares_on_grid(n):
    """
    Calculates the total number of squares on an n x n grid.

    The number of squares is given by the formula: Sum_{m=1 to n} m^2 * (n-m+1).
    This function demonstrates the calculation by printing each term of the sum
    and then the final result.
    """
    if not isinstance(n, int) or n < 1:
        print("Please provide a positive integer for n.")
        return

    total_squares = 0
    equation_terms = []

    # The summation is from m=1 to n
    for m in range(1, n + 1):
        # From the problem, the term is a^2 * b
        # We found that a = m and b = n - m + 1
        a = m
        b = n - m + 1
        
        term_value = (a**2) * b
        total_squares += term_value
        
        # This part fulfills the "output each number in the final equation" requirement.
        # It shows a^2 and b for each m.
        equation_terms.append(f"{a}^2 * {b}")

    # Joining the terms to form the full equation string
    equation_str = " + ".join(equation_terms)
    
    print(f"For an {n}x{n} grid, the total number of squares is calculated as:")
    print(f"{equation_str} = {total_squares}")


# You can change the value of n here to see the result for different grid sizes.
n_value = 4
solve_squares_on_grid(n_value)