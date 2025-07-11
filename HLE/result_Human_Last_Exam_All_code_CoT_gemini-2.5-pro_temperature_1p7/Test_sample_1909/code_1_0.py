import math

def calculate_expected_ratio(n_max):
    """
    Calculates the expected number of remaining numbers, E_n, and the ratio E_n/n.

    The function uses the recurrence relation:
    (n-1) * E_n = (n-2) * E_{n-1} + 2 * E_{n-2}
    with base cases E_0 = 0 and E_1 = 1.
    """
    if n_max < 1:
        print("Please provide a maximum n of at least 1.")
        return

    # Initialize a list to store the values of E_n.
    # E[n] will store the value of E_n.
    E = [0.0] * (n_max + 1)
    
    if n_max >= 1:
        E[1] = 1.0

    # Calculate E_n for n from 2 to n_max
    for n in range(2, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # Get the final calculated ratio for n_max
    calculated_ratio = E[n_max] / n_max if n_max > 0 else 0

    # The theoretical limit is e^(-2)
    theoretical_limit = math.exp(-2)
    
    # We are asked to output the numbers in the final equation.
    # The final equation for the limit is e^(-2).
    # The base is e and the exponent is -2.
    print(f"The equation for the limit is: e^(-2)")
    print(f"The base 'e' is approximately: {math.e}")
    print(f"The exponent is: -2")
    print("-" * 20)
    print(f"Numerical result for n = {n_max}:")
    print(f"E_{n_max} / {n_max} = {calculated_ratio:.8f}")
    print(f"Theoretical limit value:")
    print(f"e^(-2) = {theoretical_limit:.8f}")

# Set a large value for n to approximate the limit
n_large = 2000
calculate_expected_ratio(n_large)