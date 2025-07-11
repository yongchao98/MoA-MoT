import math

def solve_p_n(n):
    """
    Calculates P(n) using the derived closed-form formula.

    The closed-form formula for P(n) is:
    P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Calculate the exponent of the numerator
    exponent = n * (n - 1) / 4

    # Calculate n!
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        print(f"Error: factorial calculation is out of range for n = {n}.")
        return

    # Calculate the numerator and denominator
    numerator = (2 * math.pi) ** exponent
    denominator = math.sqrt(n_factorial)
    
    # Calculate the final result
    result = numerator / denominator

    # Output the result with the formula
    print(f"For n = {n}, the formula is:")
    print(f"P({n}) = (2 * pi)^({n} * ({n} - 1) / 4) / sqrt({n}!)")
    print(f"P({n}) = (2 * pi)^({exponent}) / sqrt({n_factorial})")
    print(f"P({n}) = {numerator} / {denominator}")
    print(f"P({n}) = {result}")

# Example: Calculate P(n) for n = 10
solve_p_n(10)