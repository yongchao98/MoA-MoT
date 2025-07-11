import numpy as np

def calculate_l(n, b):
    """
    Calculates the value of l(n, b) using the derived closed-form formula.

    The formula is: l(n,b) = 2 * (3 - (n+4)*b + 2*b^2) / (1 - b^2)
    """
    if not isinstance(n, int) or n < 10:
        raise ValueError("n must be an integer greater than or equal to 10.")
    if not isinstance(b, float) or not (0 < b < 1):
        raise ValueError("b must be a float in the interval (0, 1).")

    # The numbers in the final equation
    c1 = 2
    c2 = 3
    c3 = 4
    c4 = 2
    c5 = 1
    
    # Calculate numerator and denominator
    numerator = c1 * (c2 - (n + c3) * b + c4 * b**2)
    denominator = c5 - b**2
    
    result = numerator / denominator
    
    print(f"For n = {n} and b = {b}:")
    # Print the equation with the specific values and numbers highlighted
    print(f"l(n,b) = {c1} * ({c2} - ({n}+{c3})*{b} + {c4}*({b}**2)) / ({c5} - {b}**2)")
    print(f"l({n},{b}) = {result}")

# Example usage with n=10 and b=0.5
n_example = 10
b_example = 0.5

calculate_l(n_example, b_example)
