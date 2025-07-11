def solve():
    """
    Calculates the limit of the expected value of the ratio of the number of
    remaining numbers to n as n approaches infinity.
    """
    # Set a large number for n to approximate the limit
    N = 20000

    # Initialize a list to store the expected values E_n
    E = [0.0] * (N + 1)

    # Base cases
    E[0] = 0.0
    if N >= 1:
        E[1] = 1.0

    # Keep track of the sum sum_{k=0}^{n-2} E_k
    current_sum = E[0] + E[1]

    # Calculate E_n for n from 2 to N using the recurrence relation
    # E_n = (2 / (n-1)) * sum_{k=0}^{n-2} E_k
    for n in range(2, N + 1):
        E[n] = 2.0 * current_sum / (n - 1)
        current_sum += E[n]

    # The ratio E_n / n for a large n gives an estimate of the limit
    limit_ratio = E[N] / N

    print(f"The recurrence relation for the expected number of remaining items E_n is:")
    print("E_n = (2 / (n-1)) * sum_{k=0}^{n-2}(E_k) for n >= 2, with E_0 = 0 and E_1 = 1.")
    print(f"\nBy computing this recurrence up to n = {N}, we can estimate the limit.")
    print(f"The calculated value for E[{N}] is {E[N]:.4f}.")
    print(f"The ratio E[{N}] / {N} is {limit_ratio:.8f}.")
    print("\nThus, the limit of the expected value of the ratio is approximately:")
    print(limit_ratio)

solve()