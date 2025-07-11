def explain_ccsd_matrix_elements():
    """
    Explains which matrix elements of the form <Phi_K|H_bar|Phi> are zero in CCSD
    and why, based on the properties of the electronic Hamiltonian.
    """

    print("Step 1: Understanding the CCSD Equations")
    print("-----------------------------------------")
    print("In the Coupled Cluster Singles and Doubles (CCSD) method, the core equations to be solved for the amplitudes are:")
    print("  <Phi_S | H_bar | Phi> = 0")
    print("  <Phi_D | H_bar | Phi> = 0")
    print("where |Phi_S> and |Phi_D> are singly and doubly excited Slater determinants, respectively,")
    print("and H_bar is the similarity-transformed Hamiltonian, H_bar = exp(-T) * H * exp(T).\n")

    print("Step 2: The Role of the Hamiltonian")
    print("------------------------------------")
    print("The crucial fact is that the electronic Hamiltonian, H, contains at most two-body terms (i.e., it accounts for the kinetic energy of electrons and the Coulombic repulsion between pairs of electrons).")
    print("This property fundamentally limits the complexity of H_bar.\n")

    print("Step 3: The Structure of the Transformed Hamiltonian (H_bar)")
    print("---------------------------------------------------------------")
    print("H_bar can be expanded using the Baker-Campbell-Hausdorff (BCH) formula. Due to H being a two-body operator, this expansion terminates exactly:")
    print("H_bar = H + [H, T] + (1/2!)[[H, T], T] + (1/3!)[[[H, T], T], T] + (1/4!)[[[[H, T], T], T], T]\n")
    print("Each term in this expansion is an effective many-body operator. A detailed analysis shows that the most complex term in H_bar is an effective four-body operator.\n")

    print("Step 4: Maximum Excitation Level")
    print("---------------------------------")
    print("When H_bar acts on the reference determinant |Phi>, it creates a state vector that is a linear combination of excited determinants.")
    print("Since H_bar contains operators of up to 4-body rank, it can create excitations of at most four electrons from the reference state.")
    print("Therefore, the state H_bar|Phi> has components from the reference up to quadruply excited determinants, but no higher:")
    print("  H_bar|Phi> = c0*|Phi> + cS*|Phi_S> + cD*|Phi_D> + cT*|Phi_T> + cQ*|Phi_Q>\n")

    print("Step 5: Conclusion")
    print("-------------------")
    print("The matrix element <Phi_K|H_bar|Phi> represents the projection of the state H_bar|Phi> onto the excited determinant |Phi_K>.")
    print("Because H_bar|Phi> has no components beyond quadruple excitations, the projection onto any higher determinant is identically zero.")
    print("\nThus, in addition to the single and double excitations that are zeroed by solving the CCSD equations,")
    print("the matrix elements <Phi_K | H_bar | Phi> are guaranteed to be zero for all K-tuply excited Slater determinants where K > 4.")
    print("\nThe specific excitation levels K are:")
    final_answer = [5, 6, 7]
    print(f"For K = {final_answer[0]} (quintuples), {final_answer[1]} (sextuples), {final_answer[2]} (septuples), and all higher excitations.")


if __name__ == '__main__':
    explain_ccsd_matrix_elements()
    final_answer_string = "quintuply, sextuply, and all higher excited Slater determinants"
    # The final answer format is specified at the end.
    print(f"\n<<<The matrix elements are zero for {final_answer_string}.>>>")
