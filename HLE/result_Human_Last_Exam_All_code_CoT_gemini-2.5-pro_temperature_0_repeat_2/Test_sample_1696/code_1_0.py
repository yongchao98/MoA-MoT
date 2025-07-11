import math

def solve_and_print(n):
    """
    Calculates P(n) using the derived closed-form formula and prints the steps.

    The closed-form formula for P(n) is:
    P(n) = (2 * pi)^(n * (n - 1) / 4) * (n!)^(-1/2)
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Step 1: Calculate the exponent of (2*pi)
    exponent_2pi = n * (n - 1) / 4

    # Step 2: Calculate n!
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        print(f"Factorial of {n} is too large to compute.")
        return

    # Step 3: Calculate the exponent of n!
    exponent_factorial = -0.5

    # Step 4: Calculate the final result
    # To avoid potential overflow with large numbers, we can use logarithms
    # log(P(n)) = (n(n-1)/4) * log(2*pi) - 0.5 * log(n!)
    log_2pi = math.log(2 * math.pi)
    log_n_factorial = math.log(n_factorial)
    
    log_p_n = exponent_2pi * log_2pi + exponent_factorial * log_n_factorial
    
    # The final result can be very large, so we print it in scientific notation if needed.
    result = math.exp(log_p_n)

    print(f"For n = {n}:")
    print(f"The closed-form formula is P(n) = (2 * pi)^(n * (n - 1) / 4) * (n!)^(-1/2)")
    print(f"Substituting n = {n}:")
    print(f"P({n}) = (2 * pi)^({n} * ({n} - 1) / 4) * ({n}!)^(-1/2)")
    print(f"P({n}) = (2 * pi)^({exponent_2pi}) * ({n_factorial})^({exponent_factorial})")
    print(f"The calculated value for P({n}) is: {result}")

# Example usage with n=10
n_value = 10
solve_and_print(n_value)