def solve_set_theory_questions():
    """
    This function provides the answers to the three theoretical questions on set families.
    The reasoning for each answer is provided in the text explanation.
    """
    
    # (a) True or False: If F is a shifted (t+1)-intersecting family, then F^(1) is also (t+2)-intersecting.
    # Reasoning: This is True. A proof by contradiction shows that if we assume there is a pair
    # F1, F2 in F^(1) with intersection size t+1, we can use the shifting property to construct
    # two sets in F whose intersection is only t, which violates the (t+1)-intersecting property of F.
    answer_a = "True"

    # (b) Must a shifted (t+1)-intersecting family F satisfy |F^(n)| >= 3 for n >= k + t + 3?
    # Reasoning: This is No. A counterexample is the family F = {{1, 2, ..., k}}.
    # This family is shifted. It is (t+1)-intersecting (assuming k >= t+1).
    # Since n > k, the set {1, ..., k} does not contain n, so F^(n) = F.
    # Thus, |F^(n)| = 1, which is not >= 3.
    answer_b = "No"

    # (c) If F and G are shifted, cross-intersecting, and F is t-intersecting, does it follow
    # that F^(n) and G^(n) are also cross-intersecting?
    # Reasoning: This is Yes. This follows directly from the definitions.
    # Let F' be in F^(n) and G' be in G^(n). By definition, F' is in F and G' is in G.
    # Since F and G are cross-intersecting, |F' intersect G'| >= 1.
    # The other conditions (shifted, t-intersecting) are not needed for this conclusion.
    answer_c = "Yes"
    
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer_string)

solve_set_theory_questions()