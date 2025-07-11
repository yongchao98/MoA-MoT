def solve_complexity_analysis():
    """
    This function performs the complexity analysis for two regimes of N and L
    and prints the result in the specified format (a,b,c).
    """

    # A complexity of Theta(N log N) is represented by (a,b,c) = (2,2,0)
    # This is derived from sqrt(N^a * (log N)^b) = N * log N => a=2, b=2.
    complexity_N_log_N = (2, 2, 0)

    # For Regime 1 (N = 2^sqrt(L)), we found the complexity to be Theta(N log N).
    # L = (log2 N)^2. The divide-and-conquer strategy with chunk size l=log2(N)
    # yields a cost of N*L/l + l*2^l = N*(log N)^2/log N + (log N)*N = 2*N*log N.
    # This is the same order as the standard comparison sort, Theta(N log N).
    result_regime1 = complexity_N_log_N

    # For Regime 2 (N = 2^((log2 L)^2)), we found the complexity to be Theta(N log N).
    # L = 2^sqrt(log2 N). The divide-and-conquer strategy is less efficient than
    # a standard comparison sort because L is too large. The cost of the H-queries
    # (N*L/l) alone exceeds N*log N for any reasonable chunk size l.
    # Thus, the optimal strategy is the comparison-only sort.
    result_regime2 = complexity_N_log_N
    
    # Final answer formatting
    r1 = result_regime1
    r2 = result_regime2
    
    final_answer = f"({r1[0]},{r1[1]},{r1[2]}),({r2[0]},{r2[1]},{r2[2]})"
    print(final_answer)

solve_complexity_analysis()