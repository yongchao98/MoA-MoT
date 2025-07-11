def solve_voa_questions():
    """
    This function provides the answers to the theoretical questions about the V(p) vertex algebra.
    """
    
    # (a) Is V(p) simple for all p in Z_{>=1}?
    # No. For p >= 2, the algebra contains singular vectors which generate proper ideals,
    # making it non-simple.
    answer_a = "No"

    # (b) If there exists a non-trivial ideal in V(p), what condition must be met?
    # The existence of a non-trivial ideal in such a VOA is fundamentally due to the
    # presence of a singular (or null) vector. This vector generates the ideal.
    # Being a module for a subalgebra is a consequence of the definition of an ideal, not the cause.
    answer_b = "Must contain a singular vector"

    # (c) Does the simplicity of V(p) imply that it is also irreducible as an L_k(sl_2)-module?
    # No. Simplicity of a VOA means it has no proper ideals. However, it can still be
    # reducible as a module over a subalgebra (like its Virasoro subalgebra).
    # For example, a simple VOA can be a direct sum of multiple irreducible Virasoro modules.
    answer_c = "No"
    
    # Print the final answer in the specified format.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_voa_questions()