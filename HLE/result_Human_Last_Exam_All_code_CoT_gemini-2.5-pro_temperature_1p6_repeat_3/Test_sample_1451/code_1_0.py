def T(n):
    """
    Calculates the minimum number of trials T(n) needed to decide if a set
    of 2n bars has an equal number of real and fake ones.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("n must be a non-negative integer")
    if n == 0:
        return 0
    if n == 1:
        return 1
    
    # For n >= 2, the number of weighings depends on the parity of n.
    # The number of Phase 1 weighings is n.
    # The maximum number of Phase 2 weighings is for the largest possible
    # even number of balanced pairs (k_max_even).
    # k_max_even is n if n is even, and n-1 if n is odd.
    # Phase 2 requires k_max_even - 1 weighings.
    # Total T(n) = n + (k_max_even - 1).
    if n % 2 == 0:  # n is even
        # T(n) = n + (n - 1) = 2n - 1
        return 2 * n - 1
    else:  # n is odd
        # T(n) = n + ((n - 1) - 1) = 2n - 2
        return 2 * n - 2

# Calculate T(n) for the given values
t2 = T(2)
t3 = T(3)
t1234 = T(1234)
t6712 = T(6712)

# Print the results separated by a comma
print(f"{t2},{t3},{t1234},{t6712}")
