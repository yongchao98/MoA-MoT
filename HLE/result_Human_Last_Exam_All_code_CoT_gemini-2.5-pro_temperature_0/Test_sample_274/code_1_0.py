def count_squares_on_grid(n):
    """
    Calculates the total number of squares on an n x n grid and prints the breakdown.

    The formula for the number of squares on an n x n grid is Sum_{m=1 to n} m^2 * (n-m+1).
    This function demonstrates the calculation by showing the values for each term in the sum.
    
    Args:
        n (int): The size of the grid. Must be a positive integer.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: Please provide a positive integer for n.")
        return

    print(f"Calculating the total number of squares on a {n}x{n} grid.")
    print("The formula is: Sum_{m=1 to n} a^2 * b")
    print(f"For this problem, a = m and b = n - m + 1.\n")

    total_squares = 0
    # The summation is from m=1 to n
    for m in range(1, n + 1):
        # For the expression Sum a^2 * b, we substitute a=m and b=n-m+1
        a = m
        b = n - m + 1
        
        term = (a**2) * b
        total_squares += term
        
        # Outputting each number in the equation for the current term
        print(f"Term for m = {m}:")
        print(f"  a = m = {a}")
        print(f"  b = n - m + 1 = {n} - {m} + 1 = {b}")
        print(f"  Equation: a^2 * b = {a}^2 * {b} = {a*a} * {b} = {term}\n")

    print(f"The sum of all terms gives the total number of squares.")
    print(f"Total number of squares = {total_squares}")

# --- Main execution ---
# You can change this value to calculate for a different grid size.
grid_size = 4
count_squares_on_grid(grid_size)