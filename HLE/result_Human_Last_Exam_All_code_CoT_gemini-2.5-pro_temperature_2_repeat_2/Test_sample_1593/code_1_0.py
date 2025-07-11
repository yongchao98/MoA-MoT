def solve_complexity():
    """
    This function calculates and prints the complexity classes for the two regimes.
    
    Based on the analysis:
    - Regime 1 (N = 2^sqrt(L)): The optimal complexity is Θ(N log N).
      This corresponds to the tuple (a=2, b=2, c=0) because:
      sqrt(N^2 * (log N)^2 * (log log N)^0) = N log N.
    
    - Regime 2 (N = 2^((log L)^2)): The optimal complexity among algorithms
      that fundamentally leverage the string structure is Θ(N * (log N)^1.5).
      This corresponds to the tuple (a=2, b=3, c=0) because:
      sqrt(N^2 * (log N)^3 * (log log N)^0) = N * (log N)^(3/2).
    """
    
    # First tuple for N = 2^sqrt(L)
    r1_a, r1_b, r1_c = 2, 2, 0
    
    # Second tuple for N = 2^((log L)^2)
    r2_a, r2_b, r2_c = 2, 3, 0
    
    # Print the answer in the required format
    # The example format is "(2,0,0),(2,1,-1)".
    final_answer = f"({r1_a},{r1_b},{r1_c}),({r2_a},{r2_b},{r2_c})"
    print(final_answer)

solve_complexity()
# The final response should wrap this with the answer marker.
# <<< (2,2,0),(2,3,0) >>>