def solve_tropical_moduli_questions():
    """
    Solves the theoretical questions about tropical moduli spaces
    and prints the answer in the required format.
    """
    
    # Part (a): Minimum number of vertices for non-empty M_trop(g,A).
    # The moduli space M_trop(g,n) is non-empty if and only if 2g - 2 + n > 0.
    # When this condition holds, a stable graph with a single vertex can always be
    # constructed (a single vertex with g loops and n legs attached).
    # The valency of this vertex is 2g + n.
    # The condition 2g - 2 + n > 0 implies 2g + n > 2, which for integers means 2g + n >= 3.
    # Thus, the single vertex is at least 3-valent, satisfying the conditions.
    # The minimum number of vertices is therefore 1.
    answer_a = "1"

    # Part (b): Is M_trop(0,A) always a simplicial fan?
    # Yes, M_trop(0,n) is a classical object in tropical geometry. It is structured
    # as a fan, and its constituent cones, corresponding to stable trees, are
    # simplicial (they are isomorphic to positive orthants).
    answer_b = "yes"
    
    # Part (c): For g > 0, is M_trop(g,A) a tropical variety? If not, is it a polyhedral complex?
    # It is a tropical variety for g=1, but not for g >= 2. So the answer to the first
    # question, which pertains to all g>0, is "no".
    # For the second part ("if not..."), the space is indeed constructed by gluing
    # polyhedral cones, so it forms a generalized polyhedral complex.
    answer_c_1 = "no"
    answer_c_2 = "yes"
    
    # Format the final answer string
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c_1}, {answer_c_2}."
    
    print(final_answer)

solve_tropical_moduli_questions()