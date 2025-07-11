import math

def calculate_t_norm_1(n):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n for even n.

    Args:
        n (int): An even, non-negative integer.

    Returns:
        float: The 1-norm of the correlation matrix T.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        raise ValueError("n must be an even, non-negative integer.")

    # Using the derived formula: 2**(n+1) + 3 - (2**(n+2) + 4) / (3**n + 1)
    
    # Calculate each part of the formula
    term1 = 2**(n + 1)
    term2 = 3
    
    numerator = 2**(n + 2) + 4
    denominator = 3**n + 1
    
    term3 = numerator / denominator
    
    result = term1 + term2 - term3
    
    # Print the step-by-step calculation
    print(f"For n = {n}:")
    print(f"The 1-norm is given by the formula: 2^(n+1) + 3 - (2^(n+2) + 4) / (3^n + 1)")
    print(f"Step 1: Calculate the terms")
    print(f"  2^({n}+1) = {term1}")
    print(f"  2^({n}+2) + 4 = {2**(n+2)} + 4 = {numerator}")
    print(f"  3^{n} + 1 = {3**n} + 1 = {denominator}")
    
    print(f"Step 2: Substitute the values into the formula")
    print(f"  Norm = {term1} + {term2} - {numerator} / {denominator}")
    print(f"       = {term1 + term2} - {term3}")
    print(f"       = {result}")
    print("-" * 20)
    return result

# --- Main execution ---
# You can change this value to any even integer
even_n = 4

try:
    final_answer = calculate_t_norm_1(even_n)
    print(f"Final answer for n={even_n}: {final_answer}")
except ValueError as e:
    print(e)
