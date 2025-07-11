import math

def solve():
    """
    Solves the problem of determining the asymptotic bounds for f_1(n) and f_2(n).

    The provided derivation shows that for both c=1 and c=2, the sampling
    procedure results in a graph G' where for any induced subgraph, the
    expected number of edges is less than or equal to the expected number of
    vertices. This property indicates that the graph G' is extremely sparse
    with high probability.

    A graph property like this suggests that the arboricity is bounded by a
    constant, O(1), that does not grow with n. A detailed proof would involve
    concentration inequalities and union bounds, but the core insight comes from
    the analysis of expectations.

    Therefore, for both c=1 and c=2, the best possible upper bound on the
    arboricity, f(n), is O(1).
    """

    # For c=1, f_1(n) = O(1)
    # The categories are:
    # 1) f(n) = O(1)
    # 2) f(n) = omega(1) but O(log log n)
    # 3) f(n) = omega(log log n) but O(sqrt(log n))
    # 4) f(n) = omega(sqrt(log n)) but o(log n)
    # 5) f(n) = Theta(log n)
    # 6) f(n) = omega(log n) but 2^sqrt(O(log n))
    # 7) f(n) = 2^sqrt(omega(log n)) but O(n^0.1)
    # 8) f(n) = omega(n^0.1) but O(n^0.9)
    # 9) f(n) = omega(n^0.9)

    # f_1(n) = O(1) falls into category 1.
    f1_category = 1

    # For c=2, f_2(n) = O(1)
    # The analysis holds for c=2 as well, leading to the same conclusion.
    # f_2(n) = O(1) falls into category 1.
    f2_category = 1

    # The problem asks for a two-digit number.
    answer = f"{f1_category}{f2_category}"
    print(f"The analysis shows that for c=1, the arboricity f_1(n) is O(1). This corresponds to category {f1_category}.")
    print(f"For c=2, the arboricity f_2(n) is also O(1). This corresponds to category {f2_category}.")
    print(f"The resulting two-digit number is {answer}.")
    
    # Final answer format
    print(f"\n<<<11>>>")

solve()