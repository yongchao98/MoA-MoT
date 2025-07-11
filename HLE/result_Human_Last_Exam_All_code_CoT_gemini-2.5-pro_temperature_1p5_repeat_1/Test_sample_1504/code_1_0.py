def solve_turan_type_problem():
    """
    This function provides the solution to the Turan-type extremal problem
    by determining the truth value of statements (a) and (b) and providing
    an expression for the upper bound in (c).
    """

    # Part (a): If G is not a union of K2's, then ex(n; G, K_1,t-ind) = Theta(n).
    # This is True. A graph that is not a union of K2's contains at least one
    # connected component with 3 or more vertices, and thus contains a path P3.
    # A graph with a super-linear number of edges will have a subgraph with a high minimum degree.
    # By Ramsey's theorem, the neighborhood of a high-degree vertex must contain a large clique
    # (since it cannot contain a large independent set, which would form an induced K_1,t).
    # This large clique structure would contain G, a contradiction. Thus, the number of edges is O(n).
    # A construction of disjoint edges shows the bound is Omega(n).
    answer_a = "True"

    # Part (b): If G ~ sK2, ex(n; sK2, K_1,t-ind) is independent of n.
    # This is True. A graph H that is sK2-free has a matching number nu(H) <= s-1.
    # A graph's vertex cover number tau(H) is bounded by 2*nu(H).
    # Thus, H has a vertex cover C with |C| <= 2*(s-1). The rest of the vertices, I = V\C,
    # form an independent set. Any edge in H has at least one endpoint in C.
    # For a vertex v in C, its neighborhood N(v) cannot have an independent set of size t.
    # Since N(v) intersected with I is an independent set, |N(v) intersect I| < t.
    # The total number of edges is e(H[C]) + e(H[C,I]), which is bounded by
    # binom(|C|, 2) + |C|*(t-1). This bound depends only on s and t, not on n.
    answer_b = "True"

    # Part (c): Upper bound for ex(n; sK2, K_1,t-ind) in terms of s and t.
    # Following the logic for (b), e(H) <= binom(2s-2, 2) + (2s-2)(t-1).
    # We will provide the simplified algebraic expression.
    # binom(2s-2, 2) = (2s-2)(2s-3)/2 = (s-1)(2s-3)
    # (2s-2)(t-1) = 2(s-1)(t-1)
    # Sum = (s-1)(2s-3) + 2(s-1)(t-1) = (s-1) * (2s - 3 + 2t - 2) = (s-1)(2s + 2t - 5).
    # To satisfy the instruction "output each number in the final equation", we will
    # format the string to be explicit about the formula.
    s = "s"
    t = "t"
    expression_c = f"({s}-1)(2*{s} + 2*{t} - 5)"
    
    # The final answer in the required format
    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{expression_c}]")

solve_turan_type_problem()