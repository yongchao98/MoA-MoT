import math

def solve():
    """
    This function implements the logic derived above to determine the query complexity.
    
    The core logic is as follows:
    1.  The problem of sorting N bitstrings of length L can be approached in two main ways:
        a) A general comparison-based sort (e.g., Mergesort) using the 'C' operation on full strings.
           This has a query complexity of Theta(N * log N).
        b) A Radix Sort, using a combination of 'H' and 'C' operations on chunks of the strings.
           The optimal complexity for this method is Theta(L * N / log N).

    2.  The overall query complexity is the minimum of these two approaches. The choice depends on
        the relationship between L and N.
        - Radix Sort is better if L < (log N)^2.
        - Comparison Sort is better if L > (log N)^2.
        - They are equivalent if L is proportional to (log N)^2.

    3.  Analyze the two given regimes:
        - Regime 1: N = 2^sqrt(L) => L = (log2 N)^2.
          This is the crossover point. The complexity is Theta(N * log N).
        - Regime 2: N = 2^((log2 L)^2) => L = 2^sqrt(log2 N).
          For large N, L > (log2 N)^2. Thus, Comparison Sort is better,
          and the complexity is Theta(N * log N).

    4.  Both regimes result in a complexity of Theta(N * log N). We convert this to the (a,b,c) notation.
        The notation represents sqrt(N^a * (log N)^b * (log log N)^c).
        Theta(N * log N) corresponds to a=2, b=2, c=0.
    """

    # For both regimes, the resulting complexity is Theta(N log N).
    # We convert this to the (a, b, c) format.
    # Complexity Q = N^1 * (log N)^1 * (log log N)^0
    # Required form Q = sqrt(N^a * (log N)^b * (log log N)^c)
    #               = N^(a/2) * (log N)^(b/2) * (log log N)^(c/2)
    # Equating exponents:
    # 1 = a/2  => a = 2
    # 1 = b/2  => b = 2
    # 0 = c/2  => c = 0
    # So the tuple for both cases is (2, 2, 0).
    
    # First case (N = 2^sqrt(L)) complexity
    a1, b1, c1 = 2, 2, 0
    
    # Second case (N = 2^((log2 L)^2)) complexity
    a2, b2, c2 = 2, 2, 0
    
    # The output format is a single string "(a1,b1,c1),(a2,b2,c2)"
    result = f"({a1},{b1},{c1}),({a2},{b2},{c2})"
    print(result)

solve()