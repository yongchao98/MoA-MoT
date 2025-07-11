def solve_squares_in_grid():
    """
    Calculates the total number of squares on an n x n grid using the formula
    Sum_{m=1 to n} a^2 * b, where a=m and b=n-m+1.

    The script prompts the user for n, then prints the breakdown of the
    summation and the final result.
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

    print(f"\nFor an {n}x{n} grid, the total number of squares is calculated by the sum:")
    print(f"Sum_{{m=1 to {n}}} (m^2 * (n-m+1))")
    print("-" * 30)

    total_squares = 0
    terms = []
    for m in range(1, n + 1):
        a = m
        b = n - m + 1
        term_value = (a**2) * b
        total_squares += term_value
        terms.append(str(term_value))
        
        print(f"m={m}: a={a}, b={b}")
        print(f"     Term = {a}^2 * {b} = {a*a} * {b} = {term_value}")

    print("-" * 30)
    if n > 0:
        equation_str = " + ".join(terms)
        print(f"Total = {equation_str} = {total_squares}")
    else:
        print(f"Total = {total_squares}")

    # Verification using the closed-form formula: n(n+1)^2(n+2)/12
    verification = n * (n + 1)**2 * (n + 2) // 12
    print(f"\nVerification with closed-form formula: {verification}")
    if total_squares == verification:
        print("The result is correct.")
    else:
        print("There was an error in the calculation.")

if __name__ == '__main__':
    solve_squares_in_grid()