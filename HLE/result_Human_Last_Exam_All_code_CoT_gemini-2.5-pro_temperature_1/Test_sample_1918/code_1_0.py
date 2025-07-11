def max_complementary_choi_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H_1.
        n (int): Dimension of the output Hilbert space H_2.
        r (int): Rank of the Choi matrix of the channel Lambda.

    Returns:
        int: The maximal rank of the Choi matrix of the complementary channel Lambda^c.
    """
    # The maximal rank is given by the formula min(n, d + r, d * r).
    # We check for valid inputs first.
    if not all(isinstance(i, int) and i > 0 for i in [d, n, r]):
        raise ValueError("Dimensions and rank must be positive integers.")
    # The rank 'r' cannot exceed the dimension of the space the Choi matrix lives in, d*n.
    if r > d * n:
        raise ValueError(f"Input rank r={r} cannot be larger than d*n={d*n}.")

    # Calculate the terms in the formula
    d_plus_r = d + r
    d_times_r = d * r

    # The maximal rank is the minimum of these three values
    max_rank = min(n, d_plus_r, d_times_r)

    # Output the explanation of the calculation
    print(f"Given parameters:")
    print(f"  Dimension of input space (d): {d}")
    print(f"  Dimension of output space (n): {n}")
    print(f"  Rank of the Choi matrix (r): {r}")
    print("\nThe maximal rank of the complementary channel's Choi matrix (r^c) is given by the formula:")
    print(f"  r^c_max = min(n, d + r, d * r)")
    print(f"Substituting the given values:")
    print(f"  r^c_max = min({n}, {d} + {r}, {d} * {r})")
    print(f"  r^c_max = min({n}, {d_plus_r}, {d_times_r})")
    print(f"The result is: {max_rank}")
    
    return max_rank

if __name__ == '__main__':
    # --- Example Calculation ---
    # Dimension of the input space of the channel
    d_val = 4
    # Dimension of the output space of the channel
    n_val = 10
    # Rank of the Choi matrix of the channel
    r_val = 3

    max_complementary_choi_rank(d_val, n_val, r_val)