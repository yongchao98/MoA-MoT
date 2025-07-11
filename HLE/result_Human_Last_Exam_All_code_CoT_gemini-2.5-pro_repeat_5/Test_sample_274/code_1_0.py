import math

def solve_squares_in_grid():
    """
    Calculates the total number of squares on an n x n grid.

    The formula used is Sum_{m=1 to n} (n-m+1)^2 * m.
    This code prompts the user for n, then calculates the sum,
    printing each term of the summation.
    """
    try:
        n_str = input("Enter the grid size n: ")
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter a positive integer.")
        return

    total_squares = 0
    
    # The formula for the number of squares is Sum_{m=1 to n} a^2 * b
    # where a = n - m + 1 and b = m.
    print(f"\nCalculating the number of squares for an {n}x{n} grid.")
    print("The formula is: sum_{m=1 to n} a^2 * b")
    print("-" * 30)

    for m in range(1, n + 1):
        # In the expression sum_{m=1 to n} a^2 * b, we have:
        a = n - m + 1
        b = m
        
        term_value = (a ** 2) * b
        total_squares += term_value
        
        print(f"For m = {m}:")
        print(f"  a = n - m + 1 = {n} - {m} + 1 = {a}")
        print(f"  b = m = {b}")
        print(f"  Term = a^2 * b = {a}^2 * {b} = {term_value}")

    print("-" * 30)
    print(f"The total number of squares is: {total_squares}")

if __name__ == '__main__':
    solve_squares_in_grid()
