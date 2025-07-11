def solve():
    """
    Calculates the query complexity for the two specified regimes.
    The analysis shows that for both regimes, the complexity is Theta(N log N).
    
    Regime 1: N = 2^sqrt(L)  => L = (log N)^2
    Complexity = min(NL/log(NL), N log N)
               = min(N(log N)^2/log(N(log N)^2), N log N)
               = min(N(log N)^2/(log N + 2loglogN), N log N)
               = Theta(N log N)

    Regime 2: N = 2^((log L)^2) => L = 2^sqrt(log N)
    Complexity = min(NL/log(NL), N log N)
               = min(N*2^sqrt(logN) / log(N*2^sqrt(logN)), N log N)
               = min(N*2^sqrt(logN) / (logN + sqrt(logN)), N log N)
    Since 2^sqrt(logN) grows much faster than (logN)^2, the first term is larger.
    Complexity = Theta(N log N)

    To convert Theta(N log N) to the (a,b,c) format:
    sqrt(N^a * (log N)^b * (log log N)^c) = N log N
    N^a * (log N)^b * (log log N)^c = (N log N)^2 = N^2 * (log N)^2
    This gives a=2, b=2, c=0.
    """
    
    # For regime N = 2^sqrt(L), the complexity is (a=2, b=2, c=0)
    regime1_a = 2
    regime1_b = 2
    regime1_c = 0

    # For regime N = 2^((log L)^2), the complexity is (a=2, b=2, c=0)
    regime2_a = 2
    regime2_b = 2
    regime2_c = 0
    
    # The final answer format is (a,b,c),(a,b,c)
    final_answer = f"({regime1_a},{regime1_b},{regime1_c}),({regime2_a},{regime2_b},{regime2_c})"
    
    print(final_answer)

solve()
print("<<<({0},{1},{2}),({3},{4},{5})>>>".format(2,2,0,2,2,0))