def explain_ccsd_matrix_elements():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory.
    """

    print("Step 1: The CCSD Amplitude Equations")
    print("In CCSD, the amplitudes for the T1 and T2 cluster operators are found by ensuring that")
    print("the projection of the similarity-transformed Hamiltonian, H_bar, onto the space of")
    print("singly and doubly excited determinants is zero. This is expressed as:")
    print("  <Phi_i^a | H_bar | Phi> = 0      (for all single excitations)")
    print("  <Phi_{ij}^{ab} | H_bar | Phi> = 0   (for all double excitations)")
    print("\nThe question asks for which *other* excited determinants this matrix element is zero.")
    print("-" * 60)

    print("Step 2: The Nature of the Transformed Hamiltonian (H_bar)")
    print("The problem specifies that the electronic Hamiltonian, H, contains up to two-body terms.")
    print("The CCSD cluster operator is T = T1 + T2, where T1 and T2 are one- and two-body excitation operators.")
    print("The similarity-transformed Hamiltonian is H_bar = exp(-T) * H * exp(T).")
    print("\nA key result from many-body theory is that for a two-body H and a T limited to T1 and T2,")
    print("the resulting H_bar operator contains connected terms of at most **four-body** rank.")
    print("-" * 60)
    
    print("Step 3: Action of H_bar on the Reference State")
    print("When an operator with a maximum rank of four acts on the reference determinant |Phi>,")
    print("it can only generate excitations of up to rank four.")
    print("This means the state vector H_bar|Phi> is a linear combination of determinants from the")
    print("reference itself up to and including quadruply-excited determinants:")
    print("  H_bar|Phi> = c0|Phi> + cS|Phi_S> + cD|Phi_D> + cT|Phi_T> + cQ|Phi_Q>")
    print("-" * 60)

    print("Step 4: The Consequence of Orthogonality")
    print("The excited Slater determinants |Phi_I> are all mutually orthogonal.")
    print("The matrix element <Phi_I | H_bar | Phi> is the projection of the vector H_bar|Phi> onto |Phi_I>.")
    print("If H_bar|Phi> has no component of a certain excitation level, the projection is zero.")
    print("\nSince H_bar|Phi> has no components from quintuply-excited determinants or higher, the")
    print("matrix elements for these determinants must be identically zero.")
    print("-" * 60)

    print("Final Answer:")
    print("The matrix elements <Phi_I | H_bar | Phi> are zero by construction for singly and doubly")
    print("excited determinants.")
    print("They are generally NON-ZERO for triply and quadruply excited determinants. The non-zero")
    print("projection on the triply-excited space is the basis for the famous CCSD(T) correction.")
    print("\nThe matrix elements are identically zero due to the operator structure for:")
    print(">>> Triply and quadruply excited determinants.")

if __name__ == '__main__':
    explain_ccsd_matrix_elements()