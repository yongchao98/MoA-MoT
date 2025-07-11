def solve_ccsd_matrix_elements():
    """
    Determines for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory, based on the fact
    that the Hamiltonian is a two-body operator.
    """

    # 1. Premises from the problem statement
    # The Hamiltonian (H) is a two-body operator.
    m_H = 2
    # The CCSD cluster operator (T = T1 + T2) has a maximum of two-body excitations.
    n_T = 2

    print("Step-by-step analysis based on the operator 'body':")
    print("-" * 70)
    print(f"Premise 1: The Hamiltonian (H) is a {m_H}-body operator.")
    print(f"Premise 2: The CCSD cluster operator (T) is at most a {n_T}-body operator.")
    print("Rule: The commutator of an m-body and n-body operator is at most (m+n-1)-body.")
    print("-" * 70)

    # 2. Calculate the maximum body of H_bar through the BCH expansion
    print("Tracking the maximum operator body through the Baker-Campbell-Hausdorff expansion of H_bar:")

    # Term 0: H
    max_body = m_H
    print(f"Term 0 (H): The maximum body is {max_body}.")

    # Term 1: [H,T]
    # The equation is: new_max_body = old_max_body + n_T - 1
    old_max_body = max_body
    max_body = old_max_body + n_T - 1
    print(f"Term 1 ([H,T]): The maximum body is ({old_max_body} + {n_T} - 1) = {max_body}.")

    # Term 2: [[H,T],T]
    old_max_body = max_body
    max_body = old_max_body + n_T - 1
    print(f"Term 2 ([[H,T],T]): The maximum body is ({old_max_body} + {n_T} - 1) = {max_body}.")

    # Term 3: [[[H,T],T],T]
    old_max_body = max_body
    max_body = old_max_body + n_T - 1
    print(f"Term 3 ([[[H,T],T],T]): The maximum body is ({old_max_body} + {n_T} - 1) = {max_body}.")

    # Term 4: [[[[H,T],T],T],T]
    old_max_body = max_body
    max_body = old_max_body + n_T - 1
    print(f"Term 4 ([[[[H,T],T],T],T]): The maximum body is ({old_max_body} + {n_T} - 1) = {max_body}.")
    print("(The expansion terminates here because H is a two-body operator)")
    print("-" * 70)

    # 3. Relate operator body to excitation level
    max_excitation_level = max_body
    print(f"The analysis shows that H_bar is at most a {max_body}-body operator.")
    print(f"A {max_body}-body operator acting on the reference |Phi> can create excitations up to level {max_excitation_level}.")
    print("Therefore, the state H_bar|Phi> has no components from determinants with more than "
          f"{max_excitation_level} excitations.")
    print("-" * 70)

    # 4. Final Conclusion
    print("Final Answer:")
    print("In CCSD, the matrix elements for singly and doubly excited determinants are set to zero to solve for the amplitudes.")
    print("Matrix elements for triply, quadruply, quintuply, and sextuply excited determinants are generally non-zero.")
    print("\nThe matrix elements <Phi_I | H_bar | Phi> are guaranteed to be zero for:")
    print(f"1. All Slater determinants with an excitation level of {max_excitation_level + 1} or higher.")
    print("2. All Slater determinants that do not share the same spin and spatial symmetries as the reference determinant |Phi>.")

solve_ccsd_matrix_elements()