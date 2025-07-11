def solve_vest_complexity():
    """
    Solves the theoretical questions about the VEST problem.
    """
    # (a) Is the VEST problem, parameterised by k, still #W[2]-hard if we restrict S 
    # to the identity matrix and require all Ti matrices to commute pairwise?
    answer_a = "No"

    # (b) Is VEST #W[1]-hard if Ti matrices are restricted to diagonal Z2-matrices 
    # with at most one non-zero entry on the diagonal?
    answer_b = "Yes"

    # (c) If Ti are restricted to have only one non-zero entry in each row, 
    # what is the complexity of the decision version of VEST?
    # Options: W[1]-hard, W[2]-hard, FPT, XP
    answer_c = "W[1]-hard"

    print(f"(a) {answer_a}")
    print(f"(b) {answer_b}")
    print(f"(c) {answer_c}")

solve_vest_complexity()