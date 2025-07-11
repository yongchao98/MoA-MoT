def explain_ccsd_matrix_elements():
    """
    Explains for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD.
    """
    print("### Step-by-Step Derivation ###")
    print("\n1. The Problem:")
    print("We want to determine for which excited Slater determinants |Phi_I> the CCSD matrix element <Phi_I | H_bar | Phi> is identically zero,")
    print("given that it's zero by construction for single (I=S) and double (I=D) excitations.")
    print("H_bar is the similarity-transformed Hamiltonian: H_bar = exp(-T) * H * exp(T).")

    print("\n2. Key Properties:")
    print("- The electronic Hamiltonian, H, contains at most two-body terms. We say its 'body-ness' is 2.")
    print("- The CCSD cluster operator, T = T1 + T2, contains at most two-body excitation operators. Its 'body-ness' is also 2.")

    print("\n3. The Baker-Campbell-Hausdorff (BCH) Expansion:")
    print("H_bar can be written as a series of nested commutators. Since H is a two-body operator, this expansion terminates exactly:")
    print("H_bar = H + [H,T] + (1/2)[[H,T],T] + (1/6)[[[H,T],T],T] + (1/24)[[[[H,T],T],T],T]")

    print("\n4. Determining the 'Body-ness' of H_bar:")
    print("We can find the maximum body-ness of H_bar by analyzing each term in the BCH expansion.")
    print("We use the rule for commutators: N_body([A_m, B_n]) <= m + n - 1, where m and n are the body-ness of operators A and B.")
    print("-" * 20)

    # Initial body-ness
    h_body = 2
    t_body = 2
    print(f"Body-ness of H = {h_body}")
    print(f"Body-ness of T = {t_body}")

    # 1st commutator
    comm_1_body = h_body + t_body - 1
    print("\nCalculating body-ness of the 1st commutator, [H,T]:")
    print(f"Body-ness([H,T]) <= Body-ness(H) + Body-ness(T) - 1")
    print(f"                 <= {h_body} + {t_body} - 1 = {comm_1_body}")

    # 2nd commutator
    comm_2_body = comm_1_body + t_body - 1
    print("\nCalculating body-ness of the 2nd commutator, [[H,T],T]:")
    print(f"Body-ness([[H,T],T]) <= Body-ness([H,T]) + Body-ness(T) - 1")
    print(f"                    <= {comm_1_body} + {t_body} - 1 = {comm_2_body}")

    # 3rd commutator
    comm_3_body = comm_2_body + t_body - 1
    print("\nCalculating body-ness of the 3rd commutator, [[[H,T],T],T]:")
    print(f"Body-ness([[[H,T],T],T]) <= Body-ness([[H,T],T]) + Body-ness(T) - 1")
    print(f"                         <= {comm_2_body} + {t_body} - 1 = {comm_3_body}")

    # 4th commutator
    comm_4_body = comm_3_body + t_body - 1
    print("\nCalculating body-ness of the 4th commutator, [[[[H,T],T],T],T]:")
    print(f"Body-ness([[[[H,T],T],T],T]) <= Body-ness([[[H,T],T],T]) + Body-ness(T) - 1")
    print(f"                              <= {comm_3_body} + {t_body} - 1 = {comm_4_body}")
    print("-" * 20)

    print(f"\nSince the BCH expansion terminates here, the maximum body-ness of any operator component in H_bar is {comm_4_body}.")

    print("\n5. Relation to Excitations:")
    print("An operator with a maximum body-ness of 'n' can connect the reference determinant |Phi> to determinants with an excitation level of at most 'n'.")
    print(f"Since H_bar is at most a {comm_4_body}-body operator, the state H_bar|Phi> is a linear combination of determinants with excitation levels from 0 (reference) up to {comm_4_body} (hextuples).")
    print("H_bar|Phi> = c_0|Phi> + c_S|Phi_S> + c_D|Phi_D> + ... + c_H|Phi_H>")
    print("where H denotes hextuple (6-fold) excitations.")

    print("\n6. Conclusion:")
    print("The matrix element <Phi_I | H_bar | Phi> is the projection of H_bar|Phi> onto the determinant |Phi_I>.")
    print("Due to the orthogonality of Slater determinants, this projection is zero if |Phi_I> is not present in the expansion of H_bar|Phi>.")
    print("Therefore, the matrix element is identically zero for any excitation level higher than 6.")
    print("\nBy construction in CCSD, the matrix elements are zero for singles and doubles.")
    print("The matrix elements are generally non-zero for triples, quadruples, pentuples, and hextuples.")
    print("\nThe matrix elements <Phi_I | H_bar | Phi> are guaranteed to be zero for all Slater determinants |Phi_I> that are:")
    print(">>> Heptuple (7-fold), Octuple (8-fold), and all higher-level excitations.")

if __name__ == '__main__':
    explain_ccsd_matrix_elements()