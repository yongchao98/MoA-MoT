def solve():
    """
    This function calculates the complexity parameters (a, b, c) for the two regimes.
    The complexity is of the form Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    """

    # Regime 1: N = 2^sqrt(L)  => L = (log N)^2
    # The optimal algorithm is a Radix-style sort with complexity O(N log N).
    # Q = N * log(N)
    # Q^2 = N^2 * (log N)^2
    # Matching with N^a * (log N)^b * (log log N)^c gives:
    a1 = 2
    b1 = 2
    c1 = 0
    regime1_result = (a1, b1, c1)

    # Regime 2: N = 2^((log L)^2) => log L = sqrt(log N)
    # The optimal algorithm is a Comparison-based sort (e.g. Mergesort)
    # with complexity O(N * log N * log L).
    # Q = N * log(N) * sqrt(log N) = N * (log N)^(3/2)
    # Q^2 = N^2 * (log N)^3
    # Matching with N^a * (log N)^b * (log log N)^c gives:
    a2 = 2
    b2 = 3
    c2 = 0
    regime2_result = (a2, b2, c2)
    
    # Printing the final result in the specified format.
    print(f"{regime1_result},{regime2_result}".replace(" ", ""))

solve()