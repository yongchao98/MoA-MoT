def solve_pohozaev_questions():
    """
    Solves the conceptual questions about the Pohozaev identity and related functionals.
    """
    
    # (a) If P(u, v) = 0, is (u, v) a critical point of J?
    # No, P(u,v)=0 is a necessary but not sufficient condition. So, False.
    answer_a = "False"

    # (b) Does a unique t > 0 exist such that (u_t, v_t) is on the Pohozaev manifold?
    # Yes, this is a standard result for projecting onto such manifolds under non-criticality assumptions.
    answer_b = "Yes"

    # (c) Must the minimiser of J on {P=0} satisfy phi''(1) < 0?
    # No, this condition is related to the Nehari manifold, not the Pohozaev manifold.
    answer_c = "No"

    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")

solve_pohozaev_questions()