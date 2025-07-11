def solve_ccsd_question():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD, given a two-body Hamiltonian.
    """

    print("### Derivation ###")
    print("1. In the Coupled Cluster Singles and Doubles (CCSD) method, the core equations used to find the amplitudes (t1, t2) are derived by ensuring that the similarity-transformed Hamiltonian, H_bar, has no projection onto the space of single and double excitations:")
    print("   <Phi_singles | H_bar | Phi> = 0")
    print("   <Phi_doubles | H_bar | Phi> = 0")
    print("   where H_bar = exp(-T) * H * exp(T) and T = T1 + T2.\n")

    print("2. The question is: for which other excitation levels (triples, quadruples, etc.) is this projection also zero?\n")

    print("3. The answer lies in the nature of the H_bar operator. The electronic Hamiltonian, H, contains at most two-body interactions. The Baker-Campbell-Hausdorff expansion of H_bar is:")
    print("   H_bar = H + [H, T] + (1/2)[[H, T], T] + ...\n")

    print("4. A fundamental result of many-body theory is that for a two-body Hamiltonian H and a cluster operator T = T1 + T2, the resulting H_bar operator contains terms of at most four-body rank. This means H_bar can simultaneously interact with a maximum of four electrons.\n")

    print("5. When an operator acts on the reference determinant |Phi>, its 'body-rank' determines the maximum level of excitation it can create. A k-body operator can create excitations of rank at most k.")
    print("   Since H_bar is at most a four-body operator, the state vector H_bar|Phi> can only contain components up to quadruply excited determinants.\n")

    print("6. This implies that the state H_bar|Phi> has strictly zero projection onto the space of any determinant with an excitation level of 5 or higher.\n")

    print("### Conclusion ###")
    print("The matrix elements <Phi_I | H_bar | Phi> are zero for all Slater determinants |Phi_I> that are:")
    print("- Quintuply excited (excitation level = 5)")
    print("- Hextuply excited (excitation level = 6)")
    print("- And all higher-level excitations up to the N-electron limit.\n")

    print("It is important to note that the projections onto the triply and quadruply excited determinants are generally NON-ZERO. These non-zero projections form the basis for higher-level CC methods like CCSD(T), CCSDT, and CCSDTQ.")
    print("\nThe final equation is of the form: <Phi_I(k) | H_bar | Phi> = 0, for excitation levels k >= 5.")


if __name__ == '__main__':
    solve_ccsd_question()
    print("\n<<<Quintuply excited determinants and higher.>>>")