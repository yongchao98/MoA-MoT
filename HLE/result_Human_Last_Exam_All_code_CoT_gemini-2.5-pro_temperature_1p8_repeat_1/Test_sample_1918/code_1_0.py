def calculate_max_complementary_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of a complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H1.
        n (int): Dimension of the output Hilbert space H2.
        r (int): Rank of the Choi matrix of the original channel.
    """
    # The parameters must satisfy d <= n + r for a channel with these properties to exist.
    # We assume the given parameters come from a valid channel.
    if not (isinstance(d, int) and isinstance(n, int) and isinstance(r, int) and d > 0 and n > 0 and r > 0):
        print("Error: Dimensions and rank must be positive integers.")
        return

    # The maximal rank of the complementary channel's Choi matrix
    # is given by the formula min(n, d + r).
    max_rc = min(n, d + r)

    print("Given dimensions d = {}, n = {}, and rank r = {}.".format(d, n, r))
    print("The maximal rank of the Choi matrix of the complementary channel (rc) is given by the formula:")
    print("rc_max = min(n, d + r)")
    print(f"rc_max = min({n}, {d} + {r})")
    print(f"rc_max = min({n}, {d+r})")
    print(f"rc_max = {max_rc}")

# Example usage with placeholder values.
# Dimension of input space
d = 4
# Dimension of output space
n = 5
# Rank of the Choi matrix of the channel
r = 3

calculate_max_complementary_rank(d, n, r)