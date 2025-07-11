import numpy as np

def calculate_l_k_n(n):
    """
    Calculates the value of l_k(n) based on the derived formula.
    
    The formula is:
    l_k(n) = (1/2) * ln(n+1) - k^2 * (2n-1)/n + (n-1) * ln(k)
    
    Args:
        n (int): The dimension, must be >= 3.
    
    Returns:
        float: The calculated value of l_k(n).
    """
    if not isinstance(n, int) or n < 3:
        raise ValueError("n must be an integer greater than or equal to 3.")

    # Calculate the constant k
    k = np.log(np.sqrt(2) + 1)
    
    # Calculate the three terms of the formula for l_k(n)
    term1 = 0.5 * np.log(n + 1)
    term2 = -k**2 * (2 * n - 1) / n
    term3 = (n - 1) * np.log(k)
    
    # Sum the terms to get the final result
    result = term1 + term2 + term3
    
    # Print the breakdown of the calculation as requested
    print(f"For n = {n}:")
    print(f"k = ln(sqrt(2) + 1) = {k}")
    print("\nBreakdown of the calculation for l_k(n):")
    print(f"Term 1 (from determinant): 1/2 * ln(n+1) = {term1}")
    print(f"Term 2 (from quadratic form): -k^2 * (2n-1)/n = {term2}")
    print(f"Term 3 (from Jacobian): (n-1) * ln(k) = {term3}")
    
    print("\nThe final equation is:")
    print(f"l_k({n}) = {term1} + ({term2}) + ({term3})")
    print(f"l_k({n}) = {result}")

    return result

if __name__ == '__main__':
    # Example calculation for n = 10
    n_value = 10
    calculate_l_k_n(n_value)
    print("-" * 20)
    # Example calculation for n = 3 to check against manual derivation
    n_value_3 = 3
    calculate_l_k_n(n_value_3)