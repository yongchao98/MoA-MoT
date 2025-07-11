import math

def calculate_l(n: int):
    """
    Calculates the exact value of the function l(n) for n >= 5.

    Args:
        n: An integer greater than or equal to 5.

    Returns:
        The exact value of l(n).
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("Input n must be an integer greater than or equal to 5.")

    # The formula for l(n) is derived as:
    # l(n) = (2*n^2 + 2 - 2*(2*n - 1)*sqrt(n^2 - n + 1)) / n^2
    
    # Let's break down the final equation into its components
    term1_numerator = 2 * n**2 + 2
    term2_factor = 2 * (2 * n - 1)
    sqrt_term = math.sqrt(n**2 - n + 1)
    numerator = term1_numerator - term2_factor * sqrt_term
    denominator = n**2
    
    result = numerator / denominator

    # Output the components of the final equation as requested
    print(f"For n = {n}, the exact value of l(n) is calculated from the equation:")
    print(f"l({n}) = ( (2 * {n}^2 + 2) - (2 * (2 * {n} - 1)) * sqrt({n}^2 - {n} + 1) ) / {n}^2")
    print(f"l({n}) = ( {term1_numerator} - {term2_factor} * sqrt({n**2 - n + 1}) ) / {denominator}")
    print(f"l({n}) = {result}")

# Example calculation for n = 10, as n must be >= 5
try:
    calculate_l(10)
except ValueError as e:
    print(e)
