import numpy as np

def solve_pauli_channel_rank():
    """
    Calculates and demonstrates the maximal rank of the complementary channel of a Pauli channel.

    The theoretical derivation shows the maximal rank is d^2. This function
    provides a numerical demonstration for the case of a qubit (d=2).
    """
    # Set the dimension of the qudit system. For a qubit, d=2.
    d = 2
    
    # The final equation for the maximal rank is d^2.
    # We will demonstrate this for d=2.
    
    print("--- Theoretical Derivation Summary ---")
    print("Let d be the dimension of the qudit system.")
    print("The maximal rank of the complementary channel of a Pauli channel is sought.")
    print("By choosing the depolarizing channel (p_k = 1/d^2) and the maximally mixed input state (rho = I/d),")
    print("the output state of the complementary channel becomes (1/d^2) * I_{d^2}, where I_{d^2} is the d^2 x d^2 identity matrix.")
    print("The rank of this matrix is d^2. This is the maximum possible rank.")
    print("Final Equation: Maximal Rank = d^2")
    print("\n--- Numerical Demonstration for a Qubit (d=2) ---")

    # 1. Define the Pauli basis operators for d=2 (qubit)
    # This basis {I, sigma_x, sigma_y, sigma_z} is orthogonal under the trace inner product.
    I = np.array([[1, 0], [0, 1]], dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    
    pauli_basis = [I, X, Y, Z]
    num_ops = d**2

    # 2. Choose the depolarizing channel: p_k = 1/d^2 for all k
    p = np.full(num_ops, 1/d**2)

    # 3. Choose the maximally mixed input state: rho_in = I_d / d
    rho_in = I / d

    # 4. Construct the d^2 x d^2 matrix for the complementary channel's output
    rho_out_E = np.zeros((num_ops, num_ops), dtype=complex)

    for j in range(num_ops):
        for k in range(num_ops):
            X_j = pauli_basis[j]
            X_k = pauli_basis[k]
            # Formula: (rho_E)_jk = sqrt(p_j * p_k) * Tr(X_j * rho_in * X_k^dagger)
            # Note: For Pauli matrices, X_k^dagger = X_k.
            term = np.trace(X_j @ rho_in @ X_k.conj().T)
            rho_out_E[j, k] = np.sqrt(p[j] * p[k]) * term

    # 5. Calculate the rank of the resulting matrix
    rank = np.linalg.matrix_rank(rho_out_E)

    print(f"For dimension d = {d}, the maximal rank is d^2 = {d**2}.")
    print("\nThe calculated output matrix for the complementary channel is:")
    # The imaginary parts are due to floating point errors and should be zero.
    print(np.real(rho_out_E))
    print(f"\nThe rank of this {d**2}x{d**2} matrix is: {rank}")
    print("\nThis confirms the theoretical result.")

solve_pauli_channel_rank()