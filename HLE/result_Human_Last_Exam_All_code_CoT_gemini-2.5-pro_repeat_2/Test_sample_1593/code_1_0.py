def solve():
    """
    This function calculates the query complexity for the two specified regimes.
    
    The analysis shows that for both regimes, an optimal general algorithm is to:
    1. Use N H-queries on the full strings to identify the M unique strings.
    2. Sort these M unique strings using Theta(M log M) C-queries.
    The total complexity is N + Theta(M log M).

    In both regimes, N < 2^L, so the number of unique strings M is at most N.
    The worst-case complexity is therefore Theta(N log N).

    We convert Theta(N log N) to the (a,b,c) notation.
    The format represents Theta(sqrt(N^a * (log N)^b * (log log N)^c)).
    This is equivalent to N^(a/2) * (log N)^(b/2) * (log log N)^(c/2).

    To match N^1 * (log N)^1:
    - a/2 = 1  => a = 2
    - b/2 = 1  => b = 2
    - c/2 = 0  => c = 0
    So the complexity class for both regimes is (2,2,0).
    """

    # For the regime N = 2^sqrt(L)
    # The complexity is Theta(N log N)
    a1 = 2
    b1 = 2
    c1 = 0
    
    # For the regime N = 2^((log_2 L)^2)
    # The complexity is also Theta(N log N)
    a2 = 2
    b2 = 2
    c2 = 0
    
    # The output format is "(a,b,c),(a,b,c)"
    # We construct the string and print it.
    answer_part1 = "({},{},{})".format(a1, b1, c1)
    answer_part2 = "({},{},{})".format(a2, b2, c2)
    final_answer = "{},{}".format(answer_part1, answer_part2)
    
    print(final_answer)

solve()