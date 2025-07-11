def solve():
    """
    This function determines the best possible upper bounds on the arboricity
    of the subsampled graph G' for c=1 and c=2.

    Let f1(n) be the bound for c=1 and f2(n) be the bound for c=2.
    """

    # Case c=1: p_u = 1/d_u
    # Analysis of the complete graph K_n shows that the number of vertices
    # in the sampled graph G', |V'|, is O(sqrt(log n)) with high probability.
    # Since G' is a clique on V' vertices, its arboricity is ceil(|V'|/2).
    # Thus, f1(n) = O(sqrt(log n)).
    # This corresponds to category 3.
    f1_category = 3

    # Case c=2: p_u = 1/d_u^2
    # The sampling is much more aggressive. For K_n, the probability of
    # keeping 2 or more vertices is O(1/n^2). With high probability,
    # the resulting graph has arboricity 0. Analysis on other graphs
    # suggests that it is very difficult to form any dense subgraph.
    # Thus, f2(n) = O(1).
    # This corresponds to category 1.
    f2_category = 1

    # The result is a two-digit number formed by the categories.
    result = f"{f1_category}{f2_category}"
    
    # The final equation is the concatenation of the two category numbers.
    # We print each digit of the final result.
    # Let's formulate it as an "equation" as hinted.
    print(f"f1_category = {f1_category}")
    print(f"f2_category = {f2_category}")
    print(f"Result = {result}")

solve()