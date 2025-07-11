def solve_squares_on_grid():
    """
    Calculates the number of squares with vertices on an n x n integer grid.
    The formula used is sum_{m=1 to n} (n-m+1)^2 * m.
    """
    try:
        n_str = input("Enter the size of the grid (n): ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    total_squares = 0
    
    # The expression is sum_{m=1 to n} a^2 * b
    # Here, a = n-m+1 and b = m.
    
    print(f"\nCalculating the total number of squares for an {n}x{n} grid.")
    print("The formula is: sum_{m=1 to n} a^2 * b")
    print("where a = n-m+1 and b = m")
    print("-" * 40)
    print("m\t a\t b\t Term (a^2 * b)")
    print("-" * 40)
    
    for m in range(1, n + 1):
        a = n - m + 1
        b = m
        term = (a ** 2) * b
        total_squares += term
        print(f"{m}\t {a}\t {b}\t {a}^2 * {b} = {term}")
        
    print("-" * 40)
    print(f"Total number of squares: {total_squares}")

    # Optional: verification using the closed-form formula n*(n+1)^2*(n+2)/12
    verification = n * (n + 1)**2 * (n + 2) / 12
    print(f"Verification with closed-form formula: {int(verification)}")

if __name__ == '__main__':
    solve_squares_on_grid()
