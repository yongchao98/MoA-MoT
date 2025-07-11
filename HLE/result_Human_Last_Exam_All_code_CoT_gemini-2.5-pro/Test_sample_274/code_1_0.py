def count_grid_squares(n):
    """
    Calculates the total number of squares on an n x n grid and displays the formula.

    An n x n grid is defined by integer points (i, j) where 0 <= i, j <= n.
    The total number of squares is given by the formula: Sum_{m=1 to n} m^2 * (n - m + 1).
    This function demonstrates the calculation for a given n.

    Args:
        n (int): The size of the grid. Must be a positive integer.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: Please provide a positive integer for n.")
        return

    total_squares = 0
    equation_parts = []

    # The formula is Sum_{m=1 to n} a^2 * b
    # where a = m and b = n - m + 1
    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term = (a**2) * b
        total_squares += term
        
        # Format the string for each term of the equation
        # This fulfills the requirement to "output each number in the final equation"
        equation_parts.append(f"{a}^2 * {b}")

    # Join the parts with " + " and append the final result
    final_equation_str = " + ".join(equation_parts)
    
    # We can also show the intermediate step
    intermediate_parts = []
    for m in range(1, n + 1):
        intermediate_parts.append(str((m**2) * (n - m + 1)))
    
    print(f"For an n x n grid where n = {n}:")
    print("\nThe formula to calculate the total number of squares is Sum_{m=1 to n} (a^2 * b), where a=m and b=n-m+1.")
    print("The final equation is:")
    print(f"{final_equation_str}")
    print(f"= {' + '.join(intermediate_parts)}")
    print(f"= {total_squares}")
    
    print(f"\nThus, the total number of squares on a {n}x{n} grid is {total_squares}.")

# --- Main execution ---
# You can change the value of n to see the result for a different grid size.
n_value = 5
count_grid_squares(n_value)
