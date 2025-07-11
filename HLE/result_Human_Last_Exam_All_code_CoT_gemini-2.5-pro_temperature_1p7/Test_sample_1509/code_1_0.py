def solve_combinatorics_questions():
    """
    Solves the three true/false, yes/no questions based on combinatorial analysis.
    """
    # Analysis for each question:
    # (a) True. Proved by contradiction using the shifted property. If two sets F1, F2 in F^(1)
    #     intersect in t+1 elements, we can shift an element from their intersection to 1 in F1
    #     to create a new set F'1. The intersection of F'1 and F2 will be t, which contradicts
    #     the (t+1)-intersection property of the parent family F.
    # (b) No. A counterexample is the family F = {{1, 2, ..., k}} with t=k-1 and n >= 2k+2.
    #     This family is shifted and (t+1)-intersecting. F^(n) contains just the one set {1,..,k},
    #     so |F^(n)| = 1, which is not >= 3.
    # (c) Yes. This follows directly from the definition. F^(n) is a subset of F, and G^(n) is
    #     a subset of G. The cross-intersection property holds for all pairs from F and G, so it
    #     must hold for all pairs from their respective sub-families. The other conditions are extraneous.

    answer_a = "True"
    answer_b = "No"
    answer_c = "Yes"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_combinatorics_questions()