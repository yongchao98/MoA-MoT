import math

def calculate_T_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn
    for a given even integer n.

    Args:
        n: A non-negative even integer.
    """
    if n < 0 or n % 2 != 0:
        print("Error: n must be a non-negative even integer.")
        return

    # The derived formula for the 1-norm of T for even n is:
    # 3 + 2**(n+1) - (4 + 2**(n+2)) / (3**n + 1)
    
    # Calculate each part of the formula
    term1 = 3
    term2 = 2**(n + 1)
    
    numerator = 4 + 2**(n + 2)
    denominator = 3**n + 1
    
    term3 = numerator / denominator
    
    result = term1 + term2 - term3

    # Output the equation with the calculated numbers
    print(f"For n = {n}:")
    print(f"The formula is: 3 + 2^(n+1) - (4 + 2^(n+2)) / (3^n + 1)")
    print(f"The calculation is: {term1} + {term2} - ({numerator} / {denominator})")
    print(f"Result: {result}")

if __name__ == '__main__':
    # You can change the value of n to any even integer.
    # We use n=2 as an example.
    n_even = 2
    calculate_T_norm(n_even)
