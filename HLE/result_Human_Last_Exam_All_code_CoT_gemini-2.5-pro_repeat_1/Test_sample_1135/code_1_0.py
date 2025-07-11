import sys

def explain_ccsd_matrix_elements():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_I|H_bar|Phi> is zero in CCSD, beyond singles and doubles.
    """
    print("In the Coupled-Cluster Singles and Doubles (CCSD) method, the equations are defined by setting the projections of the similarity-transformed Hamiltonian onto the single and double excitation manifolds to zero:")
    print("  <Φ_S|H_bar|Φ> = 0")
    print("  <Φ_D|H_bar|Φ> = 0")
    print("-" * 60)
    print("The question is for which *other* excited determinants this matrix element is zero.")
    print("\nThe key is that the electronic Hamiltonian (H) contains at most two-body terms. This has a profound consequence on the structure of the transformed Hamiltonian (H_bar).")
    print("\nWhen H_bar acts on the reference determinant |Φ>, it can only generate excitations up to a certain level. Specifically, the state H_bar|Φ> contains components of single, double, triple, and quadruple excitations, but nothing higher.")
    print("\nThis means that the projections of H_bar|Φ> onto determinants with an excitation level of five or more are identically zero.")
    print("\nTherefore, the matrix elements are also zero for:")

    # The instruction "output each number in the final equation" is interpreted
    # as listing the excitation levels for which the equation <Φ_I|H_bar|Φ> = 0 holds.
    excitation_levels = {
        5: "Quintuply",
        6: "Hexatuply"
    }
    for level, name in excitation_levels.items():
        print(f"  - {name} ({level}-fold) excited determinants")
    print("  - And all higher excited determinants (7-fold, 8-fold, etc.).")
    print("-" * 60)
    print("It is important to note that for triply and quadruply excited determinants, these matrix elements are generally NOT zero. These non-zero values are the source of the leading error in the CCSD approximation.")


if __name__ == '__main__':
    explain_ccsd_matrix_elements()
