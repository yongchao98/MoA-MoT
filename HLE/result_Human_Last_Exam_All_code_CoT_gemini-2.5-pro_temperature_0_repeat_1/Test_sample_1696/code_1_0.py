import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.
    P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    Args:
        n (int): A positive integer.
        
    Returns:
        float: The calculated value of P(n).
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    # Calculate the exponent of (2*pi)
    exponent = n * (n - 1) / 4.0
    
    # Calculate n!
    n_factorial = math.factorial(n)
    
    # Calculate the numerator: (2*pi)^exponent
    try:
        numerator = math.pow(2 * math.pi, exponent)
    except OverflowError:
        return float('inf') # Return infinity if the number is too large
        
    # Calculate the denominator: sqrt(n!)
    denominator = math.sqrt(n_factorial)
    
    # Print the formula with the calculated numbers
    print(f"For n = {n}:")
    print(f"The closed-form formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print(f"P({n}) = (2*pi)^({exponent}) / sqrt({n_factorial})")
    
    # Calculate the final result
    result = numerator / denominator
    
    return result

# Example usage:
n_value = 5
try:
    result = calculate_P(n_value)
    print(f"The value of P({n_value}) is: {result}")
except (ValueError, OverflowError) as e:
    print(e)
