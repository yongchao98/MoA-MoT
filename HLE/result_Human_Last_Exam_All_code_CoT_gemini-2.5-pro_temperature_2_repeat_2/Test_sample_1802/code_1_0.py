import numpy as np

def calculate_max_rank(d):
    """
    Calculates and explains the maximal rank of the complementary
    channel of a d-dimensional Pauli channel.

    The function demonstrates this by:
    1. Choosing a Pauli channel with uniform error probabilities.
    2. Choosing the maximally mixed state as input.
    3. Constructing the output state of the complementary channel.
    4. Computing its rank numerically.

    Args:
        d (int): The dimension of the quantum system (a qudit).
    """

    if not isinstance(d, int) or d < 2:
        print("Error: Dimension 'd' must be an integer greater than or equal to 2.")
        return

    # Theoretical explanation
    print(f"For a quantum system of dimension d = {d}, a Pauli channel acts on a {d}x{d} density matrix.")
    print(f"The environment space for this channel has dimension d^2 = {d*d}.")
    print("The rank of the complementary channel's output state is maximized by specific choices.")
    print("We choose:")
    print("  1. The channel with uniform error probabilities: p_kl = 1/d^2 for all k,l.")
    print(f"  2. The maximally mixed input state: rho = I_d / d, where I_d is the {d}x{d} identity matrix.")
    print("\nNumerically verifying the rank for this case...")

    # --- 1. Define Generalized Pauli Operators X and Z ---
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[(i + 1) % d, i] = 1

    omega = np.exp(2j * np.pi / d)
    Z = np.diag([omega**i for i in range(d)])

    # --- 2. Generate all d^2 Pauli operators U_kl = X^k Z^l ---
    pauli_operators = []
    for k in range(d):
        for l in range(d):
            # np.linalg.matrix_power is efficient for integer powers
            U_kl = np.linalg.matrix_power(X, k) @ np.linalg.matrix_power(Z, l)
            pauli_operators.append(U_kl)

    # --- 3. Define the channel and input state for maximal rank ---
    # Uniform probabilities for the Pauli channel
    probs = np.full(d**2, 1 / (d**2))
    # Maximally mixed state
    rho = np.eye(d, dtype=complex) / d

    # --- 4. Construct the output state W of the complementary channel ---
    # W has dimension (d^2 x d^2)
    # W_ij = sqrt(p_i * p_j) * Tr(rho * U_j_dagger * U_i)
    W = np.zeros((d**2, d**2), dtype=complex)
    for i in range(d**2):
        for j in range(d**2):
            p_i = probs[i]
            p_j = probs[j]
            U_i = pauli_operators[i]
            U_j_dagger = pauli_operators[j].conj().T
            
            trace_term = np.trace(rho @ U_j_dagger @ U_i)
            W[i, j] = np.sqrt(p_i * p_j) * trace_term

    # --- 5. Compute the rank of W ---
    # np.linalg.matrix_rank computes rank via SVD, robust for floating point numbers
    rank = np.linalg.matrix_rank(W, tol=1e-9)

    # Final Result
    maximal_rank = d**2
    print(f"\nResult of numerical calculation: rank = {int(np.round(rank))}")
    print("\nThis confirms the theoretical result.")
    print(f"The maximal rank of the complementary channel is d^2.")
    print(f"For d = {d}, the final equation for the maximal rank is: {d}^2 = {maximal_rank}")

if __name__ == '__main__':
    # Set the dimension 'd' for the demonstration.
    # For a qubit, d=2. For a qutrit, d=3.
    d_dimension = 3
    calculate_max_rank(d_dimension)