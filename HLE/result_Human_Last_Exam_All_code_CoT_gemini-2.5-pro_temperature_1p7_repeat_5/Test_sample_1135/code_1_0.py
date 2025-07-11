import sys

def solve_task():
    """
    This script determines for which excited Slater determinants the matrix element
    <Phi_I | H_bar | Phi> is zero in CCSD theory.

    The logic is based on the operator rank in the Baker-Campbell-Hausdorff (BCH)
    expansion of the similarity-transformed Hamiltonian H_bar = exp(-T) H exp(T).
    """

    # 1. Define the initial operator ranks (body-ness).
    # H is the electronic Hamiltonian, which contains up to two-body terms.
    h_rank = 2
    # In CCSD, the cluster operator T = T1 + T2, contains up to two-body excitation operators.
    t_rank = 2

    print("Step 1: Define the ranks of the fundamental operators.")
    print(f"The rank of the Hamiltonian H is {h_rank} (it's a two-body operator).")
    print(f"The maximum rank of the CCSD cluster operator T is {t_rank} (from T2).")
    print("-" * 50)

    # 2. Analyze the BCH expansion. It terminates at the 4th commutator for a 2-body H.
    # We calculate the maximum rank of each term in the expansion:
    # H_bar = H + [H,T] + (1/2!)[[H,T],T] + ...
    print("Step 2: Analyze the Baker-Campbell-Hausdorff expansion for H_bar.")
    print("The rank of a commutator [A, B] is at most rank(A) + rank(B) - 1.")

    # Calculate ranks term by term
    print("\nCalculating ranks of the nested commutators:")
    # Term 0: H
    current_max_rank = h_rank
    print(f"Rank of H: {current_max_rank}")

    # Term 1: [H, T]
    rank_c1 = h_rank + t_rank - 1
    current_max_rank = max(current_max_rank, rank_c1)
    print(f"Rank of [H, T] <= {h_rank} + {t_rank} - 1 = {rank_c1}")

    # Term 2: [[H, T], T]
    rank_c2 = rank_c1 + t_rank - 1
    current_max_rank = max(current_max_rank, rank_c2)
    print(f"Rank of [[H, T], T] <= {rank_c1} + {t_rank} - 1 = {rank_c2}")

    # Term 3: [[[H, T], T], T]
    rank_c3 = rank_c2 + t_rank - 1
    current_max_rank = max(current_max_rank, rank_c3)
    print(f"Rank of [[[H, T], T], T] <= {rank_c2} + {t_rank} - 1 = {rank_c3}")

    # Term 4: [[[[H, T], T], T], T]
    rank_c4 = rank_c3 + t_rank - 1
    current_max_rank = max(current_max_rank, rank_c4)
    print(f"Rank of [[[[H, T], T], T], T] <= {rank_c3} + {t_rank} - 1 = {rank_c4}")
    print("The BCH expansion terminates here for a two-body Hamiltonian.")
    print("-" * 50)


    # 3. Final Conclusion
    print("Step 3: Draw conclusions based on the maximum operator rank.")
    print(f"The maximum rank of any operator term in the expansion of H_bar is {current_max_rank}.")
    print("\nA k-body operator can connect the reference determinant |Phi> to determinants with at most k excitations.")
    print(f"Therefore, the CCSD transformed Hamiltonian H_bar can create up to {current_max_rank}-fold excitations from |Phi>.")
    print("\nThis implies that the matrix element <Phi_I | H_bar | Phi> can be non-zero for any excitation I")
    print(f"from triples (3) up to hextuples ({current_max_rank}).")
    
    # By definition of CCSD, the projection onto singles and doubles is zero.
    print("\nBy definition, <Phi_I | H_bar | Phi> = 0 for singly (I=1) and doubly (I=2) excited determinants.")

    print(f"\nCrucially, for any excitation level I > {current_max_rank}, the matrix element <Phi_I | H_bar | Phi> is strictly zero")
    print("because there is no operator in H_bar that can connect |Phi> to such a highly excited state.")
    print("\nConclusion:")
    print("The 'other' excited Slater determinants for which the matrix element is zero are those with an excitation level of")
    print(f"{current_max_rank + 1} or higher.")


solve_task()