def solve_complexity():
    """
    This function implements the logic to determine the complexity class (a,b,c)
    for the two given regimes.
    """
    
    # As determined by the analysis, the complexity for both regimes is Theta(N log N).
    # We need to find a, b, c such that:
    # sqrt(N^a * (log N)^b * (log log N)^c) = N * log N
    #
    # Squaring both sides gives:
    # N^a * (log N)^b * (log log N)^c = N^2 * (log N)^2
    #
    # By comparing the exponents for N, log N, and log log N, we can determine the values.
    
    # For N: N^a = N^2  => a = 2
    # For log N: (log N)^b = (log N)^2 => b = 2
    # For log log N: (log log N)^c = (log log N)^0 => c = 0
    
    a = 2
    b = 2
    c = 0
    
    # The result is the same for both regimes given in the problem.
    case1_abc = (a, b, c)
    case2_abc = (a, b, c)
    
    print("Based on the analysis, the complexity in both regimes is Theta(N log N).")
    print("This corresponds to the following (a,b,c) parameters:")
    
    print("\nFor case 1 (N = 2^sqrt(L)):")
    print(f"a = {case1_abc[0]}")
    print(f"b = {case1_abc[1]}")
    print(f"c = {case1_abc[2]}")
    
    print("\nFor case 2 (N = 2^((log L)^2)):")
    print(f"a = {case2_abc[0]}")
    print(f"b = {case2_abc[1]}")
    print(f"c = {case2_abc[2]}")
    
    final_answer_string = f"({case1_abc[0]},{case1_abc[1]},{case1_abc[2]}),({case2_abc[0]},{case2_abc[1]},{case2_abc[2]})"
    print(f"\nFinal answer in the required format: {final_answer_string}")

solve_complexity()