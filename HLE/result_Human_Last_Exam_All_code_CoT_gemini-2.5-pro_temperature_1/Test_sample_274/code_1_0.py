def solve_squares_in_grid(n):
    """
    Calculates the total number of squares on an n x n grid and
    demonstrates the formula used.

    The formula is: Sum_{m=1 to n} a^2 * b
    where a = m and b = n - m + 1.
    """
    print(f"To find the number of squares on an {n}x{n} grid, we use the formula:")
    print(f"Total = Sum_{{m=1 to {n}}} (a^2 * b), where a = m and b = n-m+1.")
    print("-" * 40)
    print(f"Calculating the terms for n = {n}:")
    print("-" * 40)

    total_squares = 0
    equation_terms = []

    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term = a**2 * b
        total_squares += term
        
        # Output the numbers for each term in the final equation
        print(f"For m = {m}:")
        print(f"  a = {a}")
        print(f"  b = {n} - {m} + 1 = {b}")
        print(f"  The term is: {a}^2 * {b} = {term}")
        print()
        
        equation_terms.append(str(term))

    print("-" * 40)
    # Output the final equation and the total
    final_equation = " + ".join(equation_terms)
    print(f"The sum of the terms is: {final_equation}")
    print(f"The total number of squares on a {n}x{n} grid is: {total_squares}")

# You can change this value to calculate for a different grid size
grid_size = 4
solve_squares_in_grid(grid_size)