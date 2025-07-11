def solve_voa_properties():
    """
    This function determines and prints the answers to the questions about the vertex algebra V(p).
    """

    # (a) Is V(p) simple for all p in Z_>=1?
    # No. For p >= 2, V(p) is not simple because it contains singular vectors that generate proper ideals.
    answer_a = "No"

    # (b) If there exists a non-trivial ideal in V(p), what condition must be met?
    # A non-trivial ideal must be an L_k(sl_2)-module and must contain a singular vector.
    # Therefore, both conditions are required.
    answer_b = "both"

    # (c) Does the simplicity of V(p) imply that it is also irreducible as an L_k(sl_2)-module?
    # Yes. For this class of algebras, any L_k(sl_2)-submodule is also an ideal.
    # If there are no non-trivial ideals (simplicity), there can be no non-trivial submodules (irreducibility).
    answer_c = "Yes"

    # Format the final answers as requested.
    final_answer_string = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print("The answers to the questions are:")
    print(final_answer_string)
    
    # As requested, output each number from the equation k = -2 + 1/p.
    # The numbers are -2 and 1.
    print("\nNumbers from the equation k = -2 + 1/p:")
    print(-2)
    print(1)

# Execute the function to get the solution.
solve_voa_properties()