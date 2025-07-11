def solve_query_complexity():
    """
    This function calculates the query complexity for two different regimes
    of N (number of bitstrings) and L (length of bitstrings).

    The complexity is expressed in the format (a,b,c), which represents
    the class Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    Method:
    1.  A general "chunking" algorithm's complexity is formulated as Q(k) = O(N*L/k + k*2^k),
        where k is the chunk size.
    2.  This is compared against a baseline O(N log N) sorting algorithm.
    3.  The complexity for each regime is determined by choosing the better of the two algorithms.

    Regime 1: N = 2^sqrt(L)  =>  L = (log N)^2
        - The chunking algorithm with an optimal k ~ log N gives a complexity of O(N log N).
        - This is the same as the baseline.
        - O(N log N) corresponds to (a=2, b=2, c=0).

    Regime 2: N = 2^((log L)^2)  =>  L = 2^sqrt(log N)
        - The chunking algorithm is less efficient than the baseline O(N log N) sort
          because L is very large.
        - Therefore, the complexity is determined by the baseline sort, O(N log N).
        - O(N log N) corresponds to (a=2, b=2, c=0).
    """

    # For Regime 1, the complexity is O(N log N)
    # Q^2 = O(N^2 * (log N)^2)
    # This corresponds to a=2, b=2, c=0
    regime1_abc = (2, 2, 0)

    # For Regime 2, the complexity is also O(N log N)
    # Q^2 = O(N^2 * (log N)^2)
    # This corresponds to a=2, b=2, c=0
    regime2_abc = (2, 2, 0)

    # Format the final answer as requested
    final_answer = f"({regime1_abc[0]},{regime1_abc[1]},{regime1_abc[2]}),({regime2_abc[0]},{regime2_abc[1]},{regime2_abc[2]})"
    print(final_answer)

solve_query_complexity()