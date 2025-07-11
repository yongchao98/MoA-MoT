def solve():
    """
    Analyzes the graph process and determines the complexity class for each scenario.
    
    The reasoning is as follows:
    1.  The life loss of a vertex `u` in a step is `Loss(u) = sum_{v~u} 1/max(d_u, d_v)`.
    2.  A vertex survives a step if its loss is less than its current life.
    3.  Long-surviving structures require a hierarchy of degrees, where low-degree nodes are protected by high-degree nodes, which are in turn protected by even higher-degree nodes.
    4.  The number of steps in such a process is related to the number of possible "layers" in the degree hierarchy, which is limited by the maximum degree `Δ` in the graph.
    5.  This leads to the hypothesis that the maximum number of steps `T` is bounded by a function of `Δ`, likely `T = O(log Δ)`.

    Applying this hypothesis to the three cases:
    -   Case 1: Forest with maximum degree `Δ <= sqrt(log n)`.
        The maximum number of steps is bounded by `T = O(log(Δ_max)) = O(log(sqrt(log n))) = O(log log n)`.
        This corresponds to category 4.
    -   Case 2: Forest with maximum degree `Δ <= log n`.
        The maximum number of steps is bounded by `T = O(log(Δ_max)) = O(log(log n))`.
        This also corresponds to category 4.
    -   Case 3: Any forest.
        The maximum degree `Δ` can be up to `n-1`. The maximum number of steps is bounded by `T = O(log(Δ_max)) = O(log(n-1)) = O(log n)`. A corresponding lower bound can be constructed, leading to `T = Theta(log n)`.
        This corresponds to category 8.

    Combining these results, the three-digit number is 448.
    """
    
    # The first digit corresponds to f1(n)
    f1_category = 4
    
    # The second digit corresponds to f2(n)
    f2_category = 4
    
    # The third digit corresponds to f3(n)
    f3_category = 8
    
    # The final answer is the concatenation of the three digits.
    final_answer = f"{f1_category}{f2_category}{f3_category}"
    
    print(final_answer)

solve()