import math

def solve_and_print(n):
    """
    Calculates P(n) using the derived closed-form formula and prints the details.
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input 'n' must be a positive integer.")
        return

    # Calculate the components of the formula
    exponent = n * (n - 1) / 4.0
    try:
        factorial_n = math.factorial(n)
        sqrt_factorial_n = math.sqrt(factorial_n)
    except ValueError:
        print(f"Error: Factorial of n={n} is too large to compute.")
        return
        
    # Calculate the final result
    # P(n) = (2*pi)^(n(n-1)/4) / sqrt(n!)
    try:
        term_2pi = (2 * math.pi)**exponent
        result = term_2pi / sqrt_factorial_n
    except OverflowError:
        print(f"Error: Calculation for n={n} resulted in a number too large to represent.")
        return

    # Print the equation with its components
    print(f"For n = {n}:")
    print(f"The closed-form formula for P(n) is (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print(f"The exponent is n*(n-1)/4 = {n}*({n}-1)/4 = {exponent}")
    print(f"The factorial is n! = {factorial_n}")
    print(f"So, P({n}) = (2*pi)^({exponent}) / sqrt({factorial_n})")
    print(f"The final calculated value is P({n}) = {result}")

# Example usage with a positive integer n. You can change this value.
# For larger n, the result can become very large.
n_input = 4
solve_and_print(n_input)