def solve_graph_theory_problem():
    """
    This function formulates and prints the solution to the user's graph theory questions.
    """
    
    # (a) True or false: If G is not a union of K2's, then ex(n;G,K_{1,t}-ind) = Theta(n).
    # This statement is True. A lower bound can be constructed from a disjoint union of cliques,
    # which has a linear number of edges and satisfies the freeness conditions.
    # The linear upper bound is a known result in extremal graph theory.
    answer_a = "True"
    
    # (b) If G ~ sK_2, is ex(n;sK_2,K_{1,t}-ind) independent of n?
    # This is True. As derived in the thinking steps, for large enough n, the extremal
    # number becomes a constant, specifically C(2s-1, 2).
    answer_b = "True"
    
    # (c) For G ~ sK_2, express the upper bound for ex(n; sK_2,K_{1,t}-ind) in terms of s and t.
    # The upper bound is C(2s-1, 2).
    # This expression can be written as (2s-1)(2s-2)/2 which simplifies to (2s-1)(s-1).
    # Expanding this gives 2s^2 - 3s + 1.
    # The prompt requests the final code output the numbers in the final equation.
    # A clear way to express this polynomial for a symbolic variable 's' is using code syntax.
    answer_c = "2*s**2 - 3*s + 1"

    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_graph_theory_problem()