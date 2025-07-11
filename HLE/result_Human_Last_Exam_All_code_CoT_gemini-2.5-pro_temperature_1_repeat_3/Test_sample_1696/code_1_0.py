import math

def calculate_P(n):
    """
    Calculates P(n) using the derived closed-form formula.

    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Calculate the components of the formula
    exponent = n * (n - 1) / 4
    n_factorial = math.factorial(n)
    
    # Print the formula with the specific numbers for the given n
    print(f"The closed-form formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print("-" * 20)
    print(f"For n = {n}:")
    
    # The instruction asks to output each number in the final equation.
    # The numbers are the exponent and the argument to the square root.
    print(f"P({n}) = (2*pi)^({exponent}) / sqrt({n_factorial})")

    # Calculate the numerical value
    # Use try-except for potential overflow with large n
    try:
        value = (2 * math.pi)**exponent / math.sqrt(n_factorial)
        print(f"Numerical value: {value}")
    except OverflowError:
        # For large n, calculate the natural logarithm instead
        log_value = exponent * math.log(2 * math.pi) - 0.5 * math.log(n_factorial)
        print(f"Numerical value is too large to represent.")
        print(f"Natural logarithm of P({n}): {log_value}")

if __name__ == '__main__':
    # Example: Calculate P(n) for a given n.
    # You can change this value to test other integers.
    n_value = 10
    calculate_P(n_value)