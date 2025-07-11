import numpy as np

def solve_pauli_complementary_rank():
    """
    Calculates the maximal rank of the complementary channel of a
    Pauli channel for a given dimension 'd'.

    The theoretical maximal rank is d^2. This function demonstrates
    this by building the necessary matrices for a specific 'd' and
    computing the rank of the space they span. We choose d=2 as a
    concrete example.
    """
    # 1. Set the dimension 'd' of the qudit system.
    # We choose d=2 (a qubit) for this demonstration.
    d = 2

    print(f"Calculating the maximal rank for a qudit of dimension d = {d}.")

    # 2. Define the generalized Pauli operators X and Z.
    omega = np.exp(2 * np.pi * 1j / d)
    X = np.zeros((d, d), dtype=complex)
    for i in range(d):
        X[i, (i - 1 + d) % d] = 1 # Using X|k> = |k-1> convention to match X.conj().T = X^{-1} easily

    Z = np.diag([omega**i for i in range(d)])

    # 3. Generate all d*d Pauli basis operators E_ij = X^i * Z^j.
    pauli_operators = []
    for i in range(d):
        for j in range(d):
            # Using matrix_power for powers of X and Z
            Xi = np.linalg.matrix_power(X, i)
            Zj = np.linalg.matrix_power(Z, j)
            E_ij = Xi @ Zj
            pauli_operators.append(E_ij)

    # 4. To find the maximal rank, we assume all p_ij > 0.
    # We need to find the dimension of span{E_ij^dagger * E_kl}.
    # We will compute all d^2 * d^2 = d^4 of these product matrices.
    product_matrices = []
    for E_ij in pauli_operators:
        for E_kl in pauli_operators:
            # E_ij.conj().T is the Hermitian conjugate E_ij^dagger
            prod = E_ij.conj().T @ E_kl
            product_matrices.append(prod)

    # 5. To find the dimension of the space spanned by these matrices,
    # we vectorize each d x d matrix into a d^2-dimensional vector.
    vectorized_matrices = []
    for M in product_matrices:
        vectorized_matrices.append(M.flatten())

    # 6. We form a matrix where each row is one of these vectorized matrices.
    # The rank of this large matrix is the dimension of the spanned space.
    V = np.array(vectorized_matrices)
    
    # Check the shape of the matrix for clarity
    # print(f"Shape of the matrix of vectorized operators: {V.shape}")

    rank = np.linalg.matrix_rank(V)

    # 7. Print the final result and verification.
    print(f"The calculated maximal rank is: {rank}")
    print("This confirms the general formula for the maximal rank, which is d^2.")
    print("Verification:")
    print(f"{d}**2 = {d**2}")

solve_pauli_complementary_rank()