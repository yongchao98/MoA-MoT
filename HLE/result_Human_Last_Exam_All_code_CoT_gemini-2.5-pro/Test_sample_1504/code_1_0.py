def solve_graph_theory_problem():
    """
    Solves the given Turan-type extremal graph theory problem.
    """
    
    # Part (a): True or false: If G is not a union of K_2's, then ex(n; G, K_{1,t}-ind) = Θ(n).
    answer_a = "True"
    
    # Part (b): True or false: If G ~ sK_2, ex(n; sK_2, K_{1,t}-ind) is independent of n.
    answer_b = "True"
    
    # Part (c): Express the upper bound for ex(n; sK_2, K_{1,t}-ind) in terms of s and t.
    # The bound is (s-1)(t-1) + C(s-1, 2).
    # C(n, k) = n! / (k! * (n-k)!)
    # C(s-1, 2) = (s-1)*(s-2)/2
    # The expression shows all the numbers in the equation as requested.
    expression_c = "(s-1)*(t-1) + (s-1)*(s-2)/2"
    
    # Format the final answer string.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {expression_c}"
    
    print(final_answer)

solve_graph_theory_problem()