import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).

    Args:
        n (int): A positive integer.

    Returns:
        float: The value of P(n).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return None

    # The components of the formula P(n) = (2*pi)^exponent / sqrt(radicand)
    
    # Calculate the exponent of (2*pi)
    exponent_numerator = n * (n - 1)
    exponent = exponent_numerator / 4
    
    # Calculate the term under the square root in the denominator
    try:
        radicand = math.factorial(n)
    except ValueError:
        print(f"Error: n={n} is too large to compute factorial.")
        return None

    print(f"For n = {n}, we calculate the closed form P(n) = (2 * pi)^(n * (n-1) / 4) / sqrt(n!)\n")
    print("The components of the formula are:")
    print(f"The exponent on (2 * pi) is: n * (n - 1) / 4 = {n} * {n-1} / 4 = {exponent}")
    print(f"The number under the square root is n! = {radicand}\n")
    
    # Calculate P(n)
    numerator = (2 * math.pi) ** exponent
    denominator = math.sqrt(radicand)
    
    result = numerator / denominator
    
    print(f"P({n}) = (2 * pi)^{exponent} / sqrt({radicand})")
    print(f"P({n}) = {numerator} / {denominator}")
    print(f"P({n}) = {result}")
    
    return result

if __name__ == '__main__':
    # Example: Calculate P(n) for n = 4
    n_value = 4
    calculate_P(n_value)
    
    print("\n" + "="*20 + "\n")

    # Example: Calculate P(n) for n = 10
    n_value = 10
    calculate_P(n_value)
