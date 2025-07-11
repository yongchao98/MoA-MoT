def solve_vest_complexity():
    """
    Solves the theoretical questions about VEST complexity.
    """
    # (a) Is the VEST problem, parameterised by k, still #W[2]-hard if we restrict S to the identity matrix
    # and require all T_i matrices to commute pairwise?
    answer_a = "No"

    # (b) Is VEST #W[1]-hard if T_i matrices are restricted to diagonal Z_2-matrices
    # with at most one non-zero entry on the diagonal?
    answer_b = "No"

    # (c) If T_i are restricted to have only one non-zero entry in each row,
    # what is the complexity of the decision version of VEST?
    # Choose from: W[1]-hard, W[2]-hard, FPT, XP
    answer_c = "W[1]-hard"

    print(f"(a) {answer_a}")
    print(f"(b) {answer_b}")
    print(f"(c) {answer_c}")

solve_vest_complexity()