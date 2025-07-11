def solve_squares_on_grid():
    """
    Calculates the total number of squares on an n x n grid and
    prints the full summation equation.
    """
    try:
        n_str = input("Enter the grid size n: ")
        n = int(n_str)
        if n <= 0:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    total_squares = 0
    sum_parts = []
    
    # The formula is Sum_{m=1 to n} (n-m+1) * m^2
    # Here a = m and b = n-m+1
    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term = b * (a ** 2)
        total_squares += term
        sum_parts.append(f"{b}*{a}^2")

    print(f"\nFor an n={n} grid, the expressions are a=m and b=n-m+1.")
    print("The summation is:")
    
    # Print the full equation
    equation_str = " + ".join(sum_parts)
    print(f"Total Squares = {equation_str}")
    
    # Print the equation with calculated numbers
    calc_parts = []
    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term = b * (a**2)
        calc_parts.append(f"{b}*{a*a}")

    final_calc_str = " + ".join(calc_parts)
    print(f"              = {final_calc_str}")

    print(f"              = {total_squares}")


# Execute the function
solve_squares_on_grid()