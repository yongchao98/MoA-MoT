def calculate_limit_ratio(n_max):
    """
    Calculates the expected number of remaining numbers E_n and the ratio E_n/n
    up to n_max using the derived recurrence relation.

    The recurrence relation is:
    E_n = (2 / (n-1)) * sum(E_k for k=0 to n-2)

    A more efficient form for calculation is used:
    (n-1)*E_n = (n-2)*E_{n-1} + 2*E_{n-2}
    """
    # Initialize a dictionary to store the computed values of E_n
    # Base cases
    E = {0: 0, 1: 1}

    # Iteratively compute E_n for n from 2 to n_max
    for n in range(2, n_max + 1):
        if n == 2:
            # For n=2, we have (1,2), must be removed. 0 survivors.
            # From recurrence: E_2 = 2/(1) * E_0 = 0
            E[n] = 0
        else:
            # Use the more efficient recurrence: (n-1)*E_n = (n-2)*E_{n-1} + 2*E_{n-2}
            # which is derived from the sum form to avoid re-summing.
            E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # Print some initial values of the ratio E_n/n
    print("First 20 values of the ratio E_n/n:")
    for n in range(1, 21):
        ratio = E[n] / n
        print(f"n={n:2d}: {ratio:.6f}")

    # The limit is the value of the ratio for a large n
    limit_approximation = E[n_max] / n_max
    print(f"\nApproximation of the limit for n = {n_max}:")
    print(f"{E[n_max]} / {n_max} = {limit_approximation}")


# Set a large value for n to get a good approximation of the limit.
# A larger number will give a better approximation but take longer to compute.
calculate_limit_ratio(n_max=10000)