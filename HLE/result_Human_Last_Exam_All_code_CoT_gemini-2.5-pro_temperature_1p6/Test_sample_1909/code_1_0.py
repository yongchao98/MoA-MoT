def solve():
    """
    Calculates the expected number of remaining numbers, E_n, and the ratio E_n/n
    for a large n to find the limit.
    """
    n_max = 2000
    # E[i] will store the expected number of remaining items for a row of size i.
    # We use a dictionary to handle non-contiguous indices if needed, and for clarity.
    E = {}
    
    # Base cases
    # E[0] is not directly used in recurrence but helpful for initialization. For an empty row, 0 numbers remain.
    E[0] = 0
    # For a row of 1, 1 number remains.
    E[1] = 1
    # For a row of 2, 0 numbers remain.
    E[2] = 0

    # Calculate E_n using the recurrence E_n = ((n-2)*E_{n-1} + 2*E_{n-2}) / (n-1)
    for n in range(3, n_max + 1):
        E[n] = ((n - 2) * E[n - 1] + 2 * E[n - 2]) / (n - 1)

    # The result is the limit of E_n / n as n -> infinity.
    # We approximate this with n = n_max.
    n = n_max
    ratio = E[n] / n
    
    print(f"For n = {n}:")
    print(f"The expected number of remaining numbers E[{n}] is approximately {E[n]:.4f}.")
    # The prompt asks to output each number in the final equation.
    print(f"The ratio is E[{n}] / {n} = {E[n]:.4f} / {n} = {ratio:.6f}")
    print("\nThis suggests the limit approaches 1/3.")

solve()