import math

def solve():
    """
    This function calculates and prints the query complexity for the two regimes.

    The query complexity for sorting N L-bit strings using the given operations is
    min(Complexity_RadixSort, Complexity_ComparisonSort).

    Complexity_RadixSort = Theta(N * L / log(N))
    Complexity_ComparisonSort = Theta(N * log(N))

    The overall complexity is Theta(N * min(L/log(N), log(N))).

    The result is expressed in the format (a,b,c), which represents the complexity class
    Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    """

    # --- Analysis for Regime 1: N = 2^sqrt(L) ---
    # This implies log2(N) = sqrt(L), so L = (log2(N))^2.
    # We compare L/log(N) with log(N).
    # L/log(N) = (log N)^2 / log N = log N.
    # Since L/log(N) = log(N), the minimum is log(N).
    # The complexity is Theta(N * log N).
    # To convert to the (a,b,c) format:
    # N * log N = sqrt(N^2 * (log N)^2)
    #           = sqrt(N^2 * (log N)^2 * (log log N)^0)
    # So, (a,b,c) = (2,2,0).
    regime1_abc = (2, 2, 0)

    # --- Analysis for Regime 2: N = 2^((log2 L)^2) ---
    # This implies log2(N) = (log2 L)^2, so log2(L) = sqrt(log2(N)).
    # This means L = 2^sqrt(log2(N)).
    # We compare L/log(N) with log(N).
    # This is equivalent to comparing L with (log N)^2.
    # We compare L = 2^sqrt(log N) with (log N)^2.
    # Let y = log N. We compare 2^sqrt(y) with y^2.
    # For large y, the exponential function 2^sqrt(y) grows much faster than the
    # polynomial function y^2.
    # Therefore, L > (log N)^2, which means L/log(N) > log(N).
    # The minimum is log(N).
    # The complexity is Theta(N * log N).
    # As before, this corresponds to (a,b,c) = (2,2,0).
    regime2_abc = (2, 2, 0)
    
    # Printing the result in the specified format "(a,b,c),(a,b,c)"
    result_string = f"({regime1_abc[0]},{regime1_abc[1]},{regime1_abc[2]}),({regime2_abc[0]},{regime2_abc[1]},{regime2_abc[2]})"
    print(result_string)

solve()