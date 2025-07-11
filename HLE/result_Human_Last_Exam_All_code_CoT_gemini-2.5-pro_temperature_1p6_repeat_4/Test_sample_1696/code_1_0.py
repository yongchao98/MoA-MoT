import math

def solve_for_p_n(n):
    """
    Calculates the value of P(n) for a positive integer n using the derived closed-form formula.
    
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    This function prints the components of the formula and the final result.
    """
    # Input validation
    if not isinstance(n, int) or n < 1:
        print("Error: Input must be a positive integer.")
        return

    # Calculate the components of the closed-form formula
    exponent = n * (n - 1) / 4
    try:
        factorial_n = math.factorial(n)
    except ValueError:
        print(f"Error: Factorial calculation for n={n} is not feasible.")
        return

    # Output the components of the formula as requested
    print(f"For n = {n}, the closed-form formula is P(n) = (2 * pi)^ (n * (n - 1) / 4) / sqrt(n!)")
    print("Substituting the value of n:")
    print(f"The exponent term n * (n - 1) / 4 is: {n} * {n - 1} / 4 = {exponent}")
    print(f"The factorial term n! is: {factorial_n}")
    print(f"The equation becomes: P({n}) = (2 * {math.pi})^({exponent}) / sqrt({factorial_n})")

    # Calculate the final result. Use logarithms for numerical stability with large n.
    try:
        # Direct calculation for smaller n
        if n > 170: # math.factorial limit might be higher, but this is safe
             raise OverflowError
        result = math.pow(2 * math.pi, exponent) / math.sqrt(factorial_n)
    except OverflowError:
        # Logarithmic calculation for larger n
        log_p_n = exponent * math.log(2 * math.pi) - 0.5 * math.lgamma(n + 1)
        result = math.exp(log_p_n)

    # Print the final answer
    print(f"The numerical result is: P({n}) = {result}")

# Set a value for n to test the function.
# You can change this value to calculate P(n) for any positive integer.
n = 5
solve_for_p_n(n)