import math

def calculate_feynman_diagrams(n):
    """
    Calculates a(n) = (2n - 1)!!, the number of non-vanishing Feynman diagrams
    of order 2n for electron or photon propagators.
    """
    if n < 1:
        print("n must be a positive integer.")
        return

    # The starting number for the double factorial
    k = 2 * n - 1
    
    # List to store the numbers in the multiplication
    factors = []
    
    # Variable to store the final product
    result = 1
    
    # Calculate the double factorial (product of all odd numbers from k down to 1)
    for i in range(k, 0, -2):
        result *= i
        factors.append(str(i))
    
    # Format the equation string
    equation = " * ".join(factors)
    
    # Print the final result and the equation
    print(f"For n = {n}, the order is 2n = {2*n}.")
    print(f"The number of diagrams a({n}) is (2*{n} - 1)!! = {k}!!")
    print("The calculation is:")
    print(f"{equation} = {result}")

# The user asked for a(3)
n_value = 3
calculate_feynman_diagrams(n_value)