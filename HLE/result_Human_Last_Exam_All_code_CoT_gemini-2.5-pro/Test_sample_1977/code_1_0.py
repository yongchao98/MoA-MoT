import math

def calculate_correlation_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state J_n with even n.

    Args:
        n: A non-negative even integer.

    Returns:
        The 1-norm value, or None if n is invalid.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: Input 'n' must be a non-negative even integer.")
        return None

    # Using integer arithmetic to maintain precision
    # Numerator: 2 * 6^n + 3^(n+1) - 2^(n+1) - 1
    term1 = 2 * (6**n)
    term2 = 3**(n + 1)
    term3 = 2**(n + 1)
    term4 = 1
    numerator = term1 + term2 - term3 - term4

    # Denominator: 1 + 3^n
    denominator = 1 + 3**n

    # The result should be an integer
    norm = numerator // denominator

    print(f"For n = {n}:")
    print(f"  Numerator   = 2 * 6^{n} + 3^{n+1} - 2^{n+1} - 1 = {term1} + {term2} - {term3} - {term4} = {numerator}")
    print(f"  Denominator = 1 + 3^{n} = {denominator}")
    print(f"  1-Norm ||T|| = Numerator / Denominator = {norm}")
    
    return norm

if __name__ == '__main__':
    # Example: Calculate the norm for n=4
    n_example = 4
    calculate_correlation_norm(n_example)
    print("-" * 20)
    # Example: Calculate the norm for n=2
    n_example_2 = 2
    calculate_correlation_norm(n_example_2)
