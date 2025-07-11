def solve_and_print_complexity():
    """
    This function calculates and prints the query complexity for sorting N bitstrings
    of length L in two different regimes, using the specified (a,b,c) notation.

    The complexity is given by Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    Methodology:
    We consider two primary algorithms:
    1.  A "Simple Sort" algorithm with complexity Q1 = Theta(N * log(N)).
        This holds when the number of unique strings is on the order of N, which is
        true for both regimes as N << 2^L. This algorithm first uses N H-queries
        to find unique strings and then N*log(N) C-queries to sort them.
    2.  A "Radix Sort" algorithm with an optimized complexity of Q2 = Theta(L*N / log(N)).
        This algorithm breaks strings into chunks and sorts based on them.

    The overall complexity is the minimum of these two, Q = min(Q1, Q2).

    Analysis for each regime:

    Regime 1: N = 2^sqrt(L)  =>  L = (log2(N))^2
    - We compare Q1 and Q2 by comparing their variable parts: log(N) vs L/log(N).
    - L/log(N) = (log2(N))^2 / log2(N) = log2(N).
    - Since log(N) and log2(N) are of the same order, Q1 and Q2 are of the same order.
    - The complexity is Theta(N*log(N)).

    Regime 2: N = 2^((log2(L))^2)  =>  L = 2^sqrt(log2(N))
    - We compare log(N) vs L/log(N), which is equivalent to comparing (log(N))^2 vs L.
    - For large N, the term L = 2^sqrt(log2(N)) grows asymptotically faster than (log2(N))^2.
    - Therefore, L/log(N) > log(N), which means Q2 > Q1.
    - The minimum complexity is Q1 = Theta(N*log(N)).

    Conversion to (a,b,c) notation:
    For a complexity of Theta(N*log(N)), we set:
    sqrt(N^a * (log N)^b * (log log N)^c) = N * log(N)
    N^a * (log N)^b * (log log N)^c = (N * log(N))^2 = N^2 * (log N)^2
    This gives a=2, b=2, c=0.

    Both regimes result in the same complexity class.
    """

    # Regime 1: Complexity is Theta(N*log(N))
    a1, b1, c1 = 2, 2, 0

    # Regime 2: Complexity is Theta(N*log(N))
    a2, b2, c2 = 2, 2, 0

    # Print the result in the specified format "(a1,b1,c1),(a2,b2,c2)"
    # The final equation is the tuple representation of the complexity class.
    # We print each number in that final representation.
    print(f"({a1},{b1},{c1}),({a2},{b2},{c2})")

solve_and_print_complexity()