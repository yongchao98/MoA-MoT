import math

def calculate_and_print_P(n):
    """
    Calculates P(n) using the derived closed-form formula and prints the details.

    The formula is P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input must be a positive integer.")
        return

    # Calculate the exponent of 2*pi
    exponent_val = n * (n - 1) / 4.0

    # Calculate n!
    # The factorial function in Python's math module handles large integers.
    n_factorial = math.factorial(n)
    
    # Calculate sqrt(n!)
    sqrt_n_factorial = math.sqrt(n_factorial)

    # Calculate the numerator (2*pi)^exponent
    base = 2 * math.pi
    numerator_val = math.pow(base, exponent_val)

    # Calculate the final value P(n)
    if sqrt_n_factorial == 0:
        print("Error: Division by zero (sqrt(n!) is zero).")
        return

    p_n_value = numerator_val / sqrt_n_factorial

    # Print the results as requested, showing each number in the equation.
    print(f"For n = {n}, the closed-form formula for P(n) is (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("\nCalculation steps:")
    # Using 'n', 'n-1', '4', '2', 'pi' as numbers in the equation
    print(f"1. The exponent is n * (n-1) / 4 = {n} * {n-1} / 4 = {exponent_val}")
    # Using 'n' as a number in the equation
    print(f"2. The factorial term is n! = {n}! = {n_factorial}")
    print(f"3. The denominator is sqrt(n!) = sqrt({n_factorial}) = {sqrt_n_factorial}")
    print(f"4. The numerator is (2 * pi)^exponent = ({base})^{exponent_val} = {numerator_val}")
    print(f"\nFinal result:")
    print(f"P({n}) = Numerator / Denominator = {numerator_val} / {sqrt_n_factorial}")
    print(f"P({n}) = {p_n_value}")

# Example usage for n = 10
n_example = 10
calculate_and_print_P(n_example)