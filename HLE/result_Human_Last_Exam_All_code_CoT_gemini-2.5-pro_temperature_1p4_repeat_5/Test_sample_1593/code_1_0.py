def solve_query_complexity():
    """
    This function calculates and prints the query complexities for the two specified regimes.
    The complexity is represented in the (a,b,c) format, which corresponds to the
    exponents in the expression: Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    """

    # For the regime N = 2^sqrt(L), the query complexity is Theta(N * (log N)^2).
    # sqrt(N^a * (log N)^b * ...) = N * (log N)^2
    # N^(a/2) * (log N)^(b/2) * ... = N^1 * (log N)^2
    # a/2 = 1  => a = 2
    # b/2 = 2  => b = 4
    # c/2 = 0  => c = 0
    regime1_a = 2
    regime1_b = 4
    regime1_c = 0

    # For the regime N = 2^((log2 L)^2), the query complexity is Theta(N * log N).
    # sqrt(N^a * (log N)^b * ...) = N * log N
    # N^(a/2) * (log N)^(b/2) * ... = N^1 * (log N)^1
    # a/2 = 1  => a = 2
    # b/2 = 1  => b = 2
    # c/2 = 0  => c = 0
    regime2_a = 2
    regime2_b = 2
    regime2_c = 0

    # The final answer should be in the format "(a1,b1,c1),(a2,b2,c2)".
    final_answer = f"({regime1_a},{regime1_b},{regime1_c}),({regime2_a},{regime2_b},{regime2_c})"
    print(final_answer)

solve_query_complexity()