import numpy as np

def get_pauli_basis(d):
    """
    Generates the generalized Pauli basis (X^a Z^b) for a d-dimensional system.
    """
    # Shift operator X
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[(i + 1) % d, i] = 1
    
    # Phase operator Z
    omega = np.exp(2j * np.pi / d)
    Z = np.diag([omega**k for k in range(d)])
    
    basis = []
    for a in range(d):
        for b in range(d):
            # Pauli operator P_ab = X^a Z^b
            P_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            basis.append(P_ab)
            
    return basis

def calculate_complementary_channel_rank(d):
    """
    Calculates the rank of the complementary channel for a random Pauli channel
    on a d-dimensional system.
    """
    # Get the Pauli basis for dimension d
    pauli_basis = get_pauli_basis(d)
    num_pauli_ops = d**2

    # Generate a random probability distribution for the Pauli channel
    # This defines a random Pauli channel
    p = np.random.rand(num_pauli_ops)
    p /= np.sum(p)

    # Construct the Choi matrix of the complementary channel
    # J = (1/d) * sum_{j,k} sqrt(p_j*p_k) * (P_k^dagger P_j)^T kron |j><k|
    # This can be constructed as a block matrix where the (j,k)-th block of size dxd
    # is (1/d) * sqrt(p_j*p_k) * (P_k^dagger P_j)^T
    
    choi_matrix_size = d * num_pauli_ops
    J_tilde = np.zeros((choi_matrix_size, choi_matrix_size), dtype=complex)

    for j in range(num_pauli_ops):
        for k in range(num_pauli_ops):
            Pj = pauli_basis[j]
            Pk = pauli_basis[k]
            
            # Calculate the d x d block matrix
            block = np.sqrt(p[j] * p[k]) * (Pk.conj().T @ Pj).T
            
            # Place the block into the larger Choi matrix
            row_start, row_end = j * d, (j + 1) * d
            col_start, col_end = k * d, (k + 1) * d
            J_tilde[row_start:row_end, col_start:col_end] = block

    J_tilde /= d

    # Calculate the rank of the Choi matrix.
    # We use a tolerance to account for potential floating point inaccuracies.
    rank = np.linalg.matrix_rank(J_tilde, tol=1e-9)
    
    return rank

# --- Main Execution ---
if __name__ == "__main__":
    # Set the dimension of the qudit
    d = 3

    print(f"For a quantum system of dimension d = {d}:")
    
    # Calculate the theoretical maximal rank
    theoretical_rank = d**2
    print("The theoretical maximal rank of the complementary channel is d^2.")
    print("The final equation for the rank is:")
    print(f"{d}^2 = {theoretical_rank}")
    print("-" * 20)
    
    # Verify with code
    print("Verifying the result for a randomly generated Pauli channel...")
    calculated_rank = calculate_complementary_channel_rank(d)
    print(f"The numerically calculated rank is: {calculated_rank}")
    
    if calculated_rank == theoretical_rank:
        print("The numerical result matches the theoretical prediction.")
    else:
        print("The numerical result does not match the theoretical prediction.")
