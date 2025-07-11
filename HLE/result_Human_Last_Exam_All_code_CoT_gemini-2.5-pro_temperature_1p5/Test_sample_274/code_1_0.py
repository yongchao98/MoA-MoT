def calculate_total_squares(n):
    """
    Calculates the total number of squares on an n x n grid.

    The formula used is Sum_{m=1 to n} m^2 * (n-m+1).
    This function prints the breakdown of the calculation for each term in the sum.

    Args:
        n: The size of the grid.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    print(f"To find the number of squares on an {n}x{n} grid, we use the expression:")
    print(f"Sum_{{m=1 to {n}}} a^2 * b")
    print("Here, a = m and b = n - m + 1.")
    print("-" * 30)
    print(f"Calculation for n = {n}:")

    total_squares = 0
    equation_terms = []
    
    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term_value = (a**2) * b
        total_squares += term_value
        
        # Output each number in the equation for each term m
        print(f"m={m}: a={a}, b={b}  =>  {a}^2 * {b} = {a**2} * {b} = {term_value}")
        equation_terms.append(str(term_value))

    print("-" * 30)
    equation_str = " + ".join(equation_terms)
    print(f"Total number of squares = {equation_str} = {total_squares}")

# --- User execution section ---
# You can change the value of n to see the result for a different grid size.
if __name__ == "__main__":
    grid_size = 4
    calculate_total_squares(grid_size)
