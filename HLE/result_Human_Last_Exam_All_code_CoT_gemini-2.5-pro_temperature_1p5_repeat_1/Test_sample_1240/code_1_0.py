def solve_lattice_questions():
    """
    This function determines the answers to the three questions based on the
    theory of d-neighbors of integer lattices. The reasoning is explained in
    the comments.
    """

    # Question 1: For n=12, can R2(M) be of type A_11?
    # Yes. The construction with d=3 and w=(1,1,...,1) yields R2(M) = A_11.
    ans1 = "Yes"

    # Question 2: For n=15, can R2(M) contain a D_7 component?
    # Yes. The construction with d=4 and w=(2,...,2,0,...,0) (seven 2s)
    # yields R2(M) = D_7 _|_ D_8, which contains a D_7 component.
    ans2 = "yes"

    # Question 3: For n=18, d=5, can R2(M) have more than one D_k component?
    # No. For d=5 (a prime > 2), having a D_k component on indices I requires
    # w_i = 0 (mod 5) for all i in I. If two such components existed on disjoint
    # index sets I_1 and I_2, all w_i for i in I_1 U I_2 would be 0 (mod 5).
    # This forces all cross-roots to also be in M, merging the components
    # into a single D_{k_1+k_2} component, a contradiction.
    ans3 = "no"

    # Format the final answer as requested.
    final_answer_string = f"(a) {ans1}; (b) {ans2}; (c) {ans3}"
    print(final_answer_string)
    
    # The final answer in the special format for the system.
    final_answer_tag = f"({ans1.lower()}; {ans2.lower()}; {ans3.lower()})"
    print(f"<<<{final_answer_tag}>>>")

solve_lattice_questions()