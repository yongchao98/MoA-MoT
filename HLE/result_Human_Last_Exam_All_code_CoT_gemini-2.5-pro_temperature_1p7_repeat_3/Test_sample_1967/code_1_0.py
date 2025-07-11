import math

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits.

    Returns:
        int: The total number of different destabilizer sets.
    """
    if n <= 0:
        print("Number of qubits (n) must be a positive integer.")
        return 0

    print(f"Calculating the number of destabilizer sets for n = {n} qubits.")
    print("-" * 50)
    print("For each of the n destabilizer generators, we have:")
    
    # 1. Phase choices
    # Each of the n destabilizers D_i can have a global phase from {1, -1, i, -i}.
    num_phase_choices = 4**n
    print(f"1. A global phase choice from {{+1, -1, +i, -i}}. This gives 4^{n} options.")
    print(f"   For n={n}, Phase choices = 4^{n} = {num_phase_choices}")
    
    # 2. Pauli string choices
    print("\n2. A Pauli string part P_i = P_i,1...P_i,n. The choices for operators are constrained:")
    
    # 2a. Diagonal Pauli choices
    # From {S_i, D_i} = {Z_i, D_i} = 0, the i-th Pauli operator of D_i must be X or Y.
    num_diagonal_choices = 2**n
    print(f"   a) The 'diagonal' Pauli operator P_i,i must anti-commute with Z.")
    print(f"      This means P_i,i must be X or Y (2 options for each D_i).")
    print(f"      For n={n}, Diagonal choices = 2^{n} = {num_diagonal_choices}")
    
    # 2b. Off-diagonal Pauli choices
    # From [S_i, D_j]=0 for i!=j, the i-th Pauli operator of D_j must commute with Z.
    # So P_j,i must be I or Z.
    # From [D_i, D_j]=0, we must have P_i,j and P_j,i be the same (both I or both Z).
    # The number of pairs (i, j) with i < j is n*(n-1)/2.
    num_off_diagonal_pairs = n * (n - 1) // 2
    num_off_diagonal_choices = 2**num_off_diagonal_pairs
    print(f"   b) The 'off-diagonal' Pauli operators P_i,j (for i!=j) must satisfy commutation relations.")
    print(f"      This leads to 2 choices (both I or both Z) for each pair of positions (i,j) and (j,i).")
    print(f"      The number of such pairs is n(n-1)/2 = {num_off_diagonal_pairs}.")
    print(f"      For n={n}, Off-diagonal choices = 2^({num_off_diagonal_pairs}) = {num_off_diagonal_choices}")
    
    # 3. Total
    total_sets = num_phase_choices * num_diagonal_choices * num_off_diagonal_choices
    
    print("-" * 50)
    print("The total number of different sets of destabilizers is the product of these choices:")
    print(f"Total = (Phase choices) x (Diagonal choices) x (Off-diagonal choices)")
    final_equation = f"{num_phase_choices} * {num_diagonal_choices} * {num_off_diagonal_choices} = {total_sets}"
    print(final_equation)
    
    return total_sets

if __name__ == '__main__':
    # You can change the number of qubits here
    num_qubits = 4
    count_destabilizer_sets(num_qubits)
