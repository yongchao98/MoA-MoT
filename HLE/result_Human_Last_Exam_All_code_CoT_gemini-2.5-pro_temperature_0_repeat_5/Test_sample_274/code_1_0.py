def solve_grid_squares(n):
    """
    Calculates the total number of squares on an n x n grid and
    shows the calculation based on the formula:
    Sum_{m=1 to n} m^2 * (n - m + 1)
    """
    print(f"To find the number of squares on an n x n grid (for n={n}), we use the expression:")
    print(f"  Sum_{{m=1 to {n}}} a^2 * b")
    print(f"where 'a' and 'b' are filled in as:")
    print(f"  a = m")
    print(f"  b = n - m + 1")
    print("\nThis gives the formula:")
    print(f"  Total Squares = Sum_{{m=1 to {n}}} (m)^2 * (n - m + 1)")
    print("-" * 40)
    print("Calculation breakdown:")

    total_squares = 0
    # Loop from m = 1 to n
    for m in range(1, n + 1):
        # In the expression sum a^2 * b, we set a=m and b=n-m+1
        a = m
        b = n - m + 1
        term = (a**2) * b
        total_squares += term
        # Print each number in the equation for the current term
        print(f"For m={m}: a={a}, b={b}  =>  ({a})^2 * ({b}) = {term}")

    print("-" * 40)
    print(f"The total number of squares for n={n} is: {total_squares}")

# We will use n=4 as an example to demonstrate the calculation.
n_value = 4
solve_grid_squares(n_value)