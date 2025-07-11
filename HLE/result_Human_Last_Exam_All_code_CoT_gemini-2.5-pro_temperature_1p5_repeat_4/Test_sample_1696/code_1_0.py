import math

def print_p_n_formula(n):
    """
    Calculates the components of the closed formula for P(n) and prints the formula.

    The closed formula for P(n) is (2*pi)^(n*(n-1)/4) * (n!)^(-1/2).
    
    Args:
        n: A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Calculate the exponent of (2*pi)
    exponent_val = n * (n - 1) / 4
    
    # Calculate n!
    n_factorial_val = math.factorial(n)
    
    # Using 'pi' as a symbolic representation in the string
    pi_symbol = 'π'
    
    # Assemble and print the final equation string with calculated numbers
    print(f"For n = {n}:")
    print(f"The formula is P(n) = (2 * {pi_symbol})^(n*(n-1)/4) / sqrt(n!)")
    print("Plugging in the value of n:")
    print(f"P({n}) = (2 * {pi_symbol})^({n}*({n}-1)/4) / sqrt({n}!)")
    print(f"P({n}) = (2 * {pi_symbol})^({exponent_val}) / sqrt({n_factorial_val})")


# Example: find the formula for n=10
n_example = 10
print_p_n_formula(n_example)
