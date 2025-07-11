def solve_limit_ratio(n_max):
    """
    Calculates the expected ratio E_n/n using the derived recurrence relation.

    The recurrence relation for the expected number of remaining items E_n is:
    n * E_{n+1} = (n-1) * E_n + 2 * E_{n-1}
    with base cases E_0 = 0 and E_1 = 1.
    """
    if n_max < 1:
        print("n must be at least 1.")
        return

    # Initialize E_0 and E_1
    e_prev = 0  # Represents E_{i-1}
    e_curr = 1  # Represents E_i

    if n_max == 1:
        print(f"E_1 / 1 = {e_curr / 1.0}")
        return

    # Iterate from n=2 up to n_max
    for n in range(2, n_max + 1):
        # Calculate E_n using the recurrence E_n = ((n-2)*E_{n-1} + 2*E_{n-2}) / (n-1)
        # This is just a rearrangement of the main recurrence
        e_next = ((n - 2) * e_curr + 2 * e_prev) / (n - 1)
        e_prev = e_curr
        e_curr = e_next
    
    ratio = e_curr / n_max
    print(f"The ratio E_n/n for n = {n_max} is approximately: {ratio}")

# Calculate for a large n to estimate the limit
solve_limit_ratio(10000)