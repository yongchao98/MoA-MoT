import sys

def explain_ccsd_zero_matrix_elements():
    """
    Explains for which excited Slater determinants <Φ_I|H_bar|Φ> is zero in CCSD theory.
    """
    print("Background: The CCSD amplitude equations")
    print("-" * 40)
    print("In the Coupled Cluster Singles and Doubles (CCSD) method, the amplitudes for the T1 and T2 cluster operators are found by ensuring that the projection of the similarity-transformed Schrödinger equation onto the space of singly and doubly excited determinants is zero.")
    print("This is expressed as:")
    print("  <Φ_S | H_bar | Φ> = 0  (for all single excitations S)")
    print("  <Φ_D | H_bar | Φ> = 0  (for all double excitations D)")
    print("where |Φ> is the reference determinant and H_bar is the similarity-transformed Hamiltonian.")
    print("\nThe question is: for which *other* excited determinants |Φ_I> is this matrix element also zero?")
    print("\n")
    
    print("Analysis: The structure of the CCSD H_bar")
    print("-" * 40)
    print("The Hamiltonian H contains at most two-body terms. The CCSD cluster operator is T = T1 + T2.")
    print("The similarity-transformed Hamiltonian is H_bar = exp(-T) * H * exp(T).")
    print("Crucially, due to the linked-cluster theorem, only terms where all operators (H, T1, T2) are diagrammatically connected contribute to the equations.")
    print("\nA key result of many-body theory is that for CCSD, this connected form of H_bar contains operators that can create excitations of at most rank four (quadruple excitations) when acting on the reference determinant |Φ>.")
    print("This means the resulting state vector H_bar|Φ> is a linear combination of states up to quadruples:")
    print("  H_bar|Φ> = c_0|Φ> + c_S|Φ_S> + c_D|Φ_D> + c_T|Φ_T> + c_Q|Φ_Q>")
    print("\n")

    print("Conclusion: Identifying the zero matrix elements")
    print("-" * 40)
    print("From the structure of H_bar|Φ>, we can determine which matrix elements are zero:")
    print(" - <Φ_S|H_bar|Φ> and <Φ_D|H_bar|Φ> are zero by definition of the CCSD equations.")
    print(" - <Φ_T|H_bar|Φ> and <Φ_Q|H_bar|Φ> are, in general, NOT zero.")
    print("\nSince H_bar|Φ> has no components from quintuple or higher excitations, the projections onto these determinants must be zero.")
    print("\nTherefore, the matrix elements <Φ_I|H_bar|Φ> are also zero for the following excited Slater determinants |Φ_I>:")
    print("  - Quintuply excited determinants (excitation level 5)")
    print("  - Hextuply excited determinants (excitation level 6)")
    print("  - ... and so on.")
    print("\nIn general, the matrix elements are zero for all Slater determinants with an excitation level of five or higher.")

if __name__ == "__main__":
    explain_ccsd_zero_matrix_elements()
    # The final answer can be summarized as:
    # All excited Slater determinants with an excitation level of five or more.
    sys.stdout.write("<<<all excited Slater determinants with an excitation level of five or higher>>>")