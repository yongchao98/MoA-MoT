def solve():
    """
    This function analyzes the query complexity for the two given regimes and prints the result.

    The analysis identifies two main algorithmic strategies:
    1.  A direct comparison-based sort (e.g., Mergesort) using operation C.
        This has a complexity of Theta(N log N) queries.
    2.  A radix-sort-style algorithm that uses operation H to group identical substrings (blocks)
        of a chosen size k, and operation C to sort the unique blocks.
        The complexity of this method is Theta(N*L/k + U*log(U)), where U is the number of unique blocks.

    For both specified regimes, the analysis shows that the optimal query complexity is Theta(N log N).

    Regime 1: N = 2^sqrt(L)  ==> L = (log N)^2
        - Comparison sort is Theta(N log N).
        - Radix sort with optimal k=log(N) is also Theta(N log N).
        - Result: Theta(N log N).

    Regime 2: N = 2^((log L)^2) ==> L = 2^sqrt(log N)
        - Comparison sort is Theta(N log N).
        - The radix sort algorithm proves to be less efficient than Theta(N log N) for any choice of k.
        - Result: Theta(N log N).

    Conversion to (a,b,c) format:
    The complexity is represented as Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    To match Theta(N log N), we can write it as Theta(sqrt((N log N)^2)) = Theta(sqrt(N^2 * (log N)^2)).
    This gives:
    a = 2
    b = 2
    c = 0
    The complexity class for both regimes is (2,2,0).
    """

    # For Regime 1, N = 2^sqrt(L), the complexity is Theta(N log N).
    # In (a,b,c) notation, this is (2,2,0).
    regime1_abc = (2, 2, 0)

    # For Regime 2, N = 2^((log L)^2), the complexity is also Theta(N log N).
    # In (a,b,c) notation, this is (2,2,0).
    regime2_abc = (2, 2, 0)

    # Format the output string as requested.
    result = f"({regime1_abc[0]},{regime1_abc[1]},{regime1_abc[2]}),({regime2_abc[0]},{regime2_abc[1]},{regime2_abc[2]})"
    print(result)

solve()