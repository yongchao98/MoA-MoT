import math

def calculate_correlation_norm(n: int):
    """
    Calculates the 1-norm of the correlation matrix T for the state Jn with an even n.
    
    Args:
        n: A non-negative even integer.
    """
    if not isinstance(n, int) or n < 0 or n % 2 != 0:
        print("Error: Input 'n' must be a non-negative even integer.")
        return

    # The problem asks to find the 1-norm of the correlation matrix T.
    # After a detailed derivation, the 1-norm ||T||_1 is given by the formula:
    # ||T||_1 = (2**(n+1)*(3**n - 1) + (3**(n+1) - 1)) / (2**(n+1)*(3**n + 1))
    
    # We are asked to output each number in the final equation.
    # We will compute the numerator and the denominator and print their components.

    # Calculate the numerator
    num_term1 = (2**(n + 1)) * (3**n - 1)
    num_term2 = 3**(n + 1) - 1
    numerator = num_term1 + num_term2

    # Calculate the denominator
    den_factor1 = 2**(n + 1)
    den_factor2 = 3**n + 1
    denominator = den_factor1 * den_factor2

    # The final result
    norm_value = numerator / denominator

    print("The symbolic formula for the 1-norm is: (2**(n+1)*(3**n - 1) + 3**(n+1) - 1) / (2**(n+1)*(3**n + 1))")
    print(f"\n--- Calculation for n = {n} ---")

    print(f"\nNumerator is the sum of two terms:")
    print(f"Term 1: 2^({n}+1) * (3^{n} - 1) = {2**(n+1)} * ({3**n} - 1) = {num_term1}")
    print(f"Term 2: 3^({n}+1) - 1 = {3**(n+1)} - 1 = {num_term2}")
    print(f"Total Numerator = {num_term1} + {num_term2} = {numerator}")

    print(f"\nDenominator is the product of two factors:")
    print(f"Factor 1: 2^({n}+1) = {den_factor1}")
    print(f"Factor 2: 3^{n} + 1 = {3**n} + 1 = {den_factor2}")
    print(f"Total Denominator = {den_factor1} * {den_factor2} = {denominator}")
    
    print(f"\nResult: ||T||_1 = {numerator} / {denominator} = {norm_value}")


# You can change the value of `n` below to any non-negative even integer.
# The problem asks for the formula for an arbitrary even n.
# We will demonstrate with n=2 and n=4 as examples.
example_n_1 = 2
calculate_correlation_norm(example_n_1)

print("\n" + "="*50 + "\n")

example_n_2 = 4
calculate_correlation_norm(example_n_2)
