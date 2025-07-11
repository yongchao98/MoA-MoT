def solve_complexity():
    """
    This function determines the query complexity for sorting bitstrings
    in two different regimes and prints the result in the specified format.

    The reasoning is as follows:
    1.  The optimal query complexity Q for this problem is found to be
        Theta(N * log(N)) in any regime where L > log2(N).
    2.  Both specified regimes, N = 2^sqrt(L) and N = 2^((log2 L)^2),
        satisfy the condition L > log2(N) for large L.
    3.  The complexity C = Theta(N * log(N)) needs to be converted to the
        (a,b,c) format, which is based on the equation:
        C^2 = N^a * (log N)^b * (log log N)^c.
    4.  Squaring our complexity gives: C^2 = (N * log N)^2 = N^2 * (log N)^2.
    5.  By comparing the exponents with the format equation, we get the values for a, b, and c.
    """

    # For a complexity of Theta(N * log(N)), we find a, b, and c.
    # C^2 = N^2 * (log N)^2 * (log log N)^0
    final_a = 2
    final_b = 2
    final_c = 0

    # Both regimes result in the same complexity and thus the same (a,b,c) tuple.
    regime1_tuple = (final_a, final_b, final_c)
    regime2_tuple = (final_a, final_b, final_c)

    # Format the final answer string as per the problem description "(a1,b1,c1),(a2,b2,c2)".
    # The instruction "output each number in the final equation" is satisfied by showing
    # the explicit values for a,b,c being used to construct the result.
    part1_str = f"({regime1_tuple[0]},{regime1_tuple[1]},{regime1_tuple[2]})"
    part2_str = f"({regime2_tuple[0]},{regime2_tuple[1]},{regime2_tuple[2]})"

    final_answer = f"{part1_str},{part2_str}"

    print(final_answer)

solve_complexity()