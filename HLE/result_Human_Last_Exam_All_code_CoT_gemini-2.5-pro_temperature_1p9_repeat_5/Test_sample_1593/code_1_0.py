def solve_query_complexity():
    """
    This function calculates and prints the query complexity for the two specified regimes.
    The derivation is explained in the text preceding this code block.
    """

    # For the regime N = 2^sqrt(L), the complexity is Theta(N log N).
    # In (a,b,c) notation, this is (2,2,0).
    regime_1_a = 2
    regime_1_b = 2
    regime_1_c = 0

    # For the regime N = 2^((log_2 L)^2), the complexity is also Theta(N log N).
    # In (a,b,c) notation, this is (2,2,0).
    regime_2_a = 2
    regime_2_b = 2
    regime_2_c = 0

    # The final answer format is (a,b,c),(a,b,c)
    print(f"({regime_1_a},{regime_1_b},{regime_1_c}),({regime_2_a},{regime_2_b},{regime_2_c})")

solve_query_complexity()