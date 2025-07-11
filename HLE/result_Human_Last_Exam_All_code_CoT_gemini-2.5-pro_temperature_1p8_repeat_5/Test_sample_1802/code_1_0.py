import numpy as np

def generalized_pauli_operators(d):
    """
    Generates the d x d generalized Pauli operators (Weyl-Heisenberg operators).
    These are defined as U_{ab} = X^a Z^b.
    """
    operators = []
    # The d-th root of unity
    omega = np.exp(2j * np.pi / d)
    
    # Shift operator X: X|j> = |j+1>
    X = np.roll(np.eye(d, dtype=complex), 1, axis=0)
    
    # Phase operator Z: Z|j> = omega^j |j>
    Z = np.diag([omega**j for j in range(d)])
    
    # Generate all d^2 operators U_{a,b} = X^a Z^b
    for a in range(d):
        for b in range(d):
            U_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            operators.append(U_ab)
            
    return operators

def main():
    """
    Calculates and demonstrates the maximal rank of the complementary channel for a Pauli channel on a qudit system.
    """
    # The dimension of the qudit system. Change this value to test for different dimensions.
    d = 3
    print(f"This script calculates the maximal rank of the complementary channel of a Pauli channel for a qudit of dimension d = {d}.")

    print("\n--- Summary of Theoretical Derivation ---")
    print("1. The rank of the complementary channel is the rank of the Gram matrix K of the channel's Kraus operators.")
    print("2. For a Pauli channel, the Gram matrix K is a diagonal matrix.")
    print("3. The diagonal entries are K_ii = d * p_i, where p_i are the probabilities defining the channel.")
    print("4. The rank of K is the number of non-zero probabilities p_i.")
    print(f"5. To maximize the rank, all d^2 = {d*d} probabilities must be non-zero.")

    # --- Numerical Verification ---
    print(f"\n--- Numerical Verification for d = {d} ---")
    
    # Number of Pauli operators
    num_ops = d**2

    # To maximize the rank, we choose all probabilities p_k to be non-zero.
    # We use the completely depolarizing channel where all p_k are equal.
    probabilities = np.full(num_ops, 1.0 / num_ops)
    print(f"We model a channel where all {num_ops} probabilities are non-zero (e.g., p_k = 1/{num_ops}).")
    
    # Step 1: Generate the generalized Pauli operators
    pauli_ops = generalized_pauli_operators(d)
    
    # Step 2: Construct the corresponding Kraus operators A_k = sqrt(p_k) * U_k
    kraus_ops = [np.sqrt(p) * U for p, U in zip(probabilities, pauli_ops)]

    # Step 3: Construct the Gram matrix K_ij = Tr(A_i^dagger A_j)
    gram_matrix = np.zeros((num_ops, num_ops), dtype=complex)
    for i in range(num_ops):
        for j in range(num_ops):
            gram_matrix[i, j] = np.trace(np.conjugate(kraus_ops[i]).T @ kraus_ops[j])

    # Step 4: Calculate the rank of the Gram matrix
    # We use a small tolerance to account for floating-point inaccuracies.
    rank = np.linalg.matrix_rank(gram_matrix, tol=1e-9)
    
    print("We now compute the rank of the corresponding Gram matrix numerically.")
    print(f"Calculated rank: {int(np.round(rank))}")
    
    # Step 5: Final conclusion and equation
    maximal_rank = d**2
    print("\nThe result confirms the theoretical derivation.")
    print("The maximal rank is given by the square of the qudit dimension d.")
    
    print("\nFinal Equation:")
    print(f"Maximal Rank = d^2")
    print(f"For d = {d}, the calculation is: Maximal Rank = {d}^2 = {maximal_rank}")

if __name__ == '__main__':
    main()
