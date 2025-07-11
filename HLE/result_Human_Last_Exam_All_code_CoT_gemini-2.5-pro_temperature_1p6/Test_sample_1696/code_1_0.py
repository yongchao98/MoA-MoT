import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.

    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    Args:
        n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    print(f"Calculating P(n) for n = {n}")
    print(f"Formula: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("-" * 30)

    # Calculate each part of the equation
    try:
        # Numerator calculation
        exponent = n * (n - 1) / 4.0
        numerator = (2 * math.pi)**exponent

        # Denominator calculation
        n_factorial = math.factorial(n)
        denominator = math.sqrt(n_factorial)

        # Final result
        result = numerator / denominator

        print(f"Exponent for the numerator: n * (n - 1) / 4 = {n} * {n-1} / 4 = {exponent}")
        print(f"Numerator: (2 * pi)^({exponent}) = {numerator}")
        print(f"Denominator: sqrt(n!) = sqrt({n_factorial}) = {denominator}")
        print("-" * 30)
        print(f"P({n}) = {numerator} / {denominator}")
        print(f"Final Result P({n}) = {result}")

    except OverflowError:
        print(f"Error: Calculation for n={n} resulted in a number too large to represent.")
    except ValueError:
        print(f"Error: factorial({n}) is not defined for non-integer or negative n.")


if __name__ == '__main__':
    # Example: Calculate P(10)
    calculate_P(10)