def solve_complexity():
    """
    Calculates and prints the query complexity for the two regimes.
    The complexity is given in the format (a,b,c), which corresponds to
    the expression Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    """

    # For the first regime, N = 2^sqrt(L).
    # The analysis shows the optimal query complexity is Theta(N * log N).
    # To match the format, we express this as:
    # N * log N = sqrt(N^2 * (log N)^2)
    # This implies a=2, b=2, c=0.
    regime1_a = 2
    regime1_b = 2
    regime1_c = 0

    # For the second regime, N = 2^((log_2 L)^2).
    # The analysis shows the optimal query complexity is Theta(N * (log N)^(3/2)).
    # To match the format, we express this as:
    # N * (log N)^(3/2) = sqrt(N^2 * (log N)^3)
    # This implies a=2, b=3, c=0.
    regime2_a = 2
    regime2_b = 3
    regime2_c = 0

    # The problem requests the final answer as a string "(a,b,c),(a,b,c)".
    # Here, we construct this string from our derived numbers.
    final_answer = f"({regime1_a},{regime1_b},{regime1_c}),({regime2_a},{regime2_b},{regime2_c})"
    
    print(final_answer)

solve_complexity()