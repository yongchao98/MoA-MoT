import numpy as np

def rank_F2(A):
    """Computes the rank of a binary matrix A over the field F_2."""
    mat = A.copy().astype(int)
    rows, cols = mat.shape
    rank = 0
    pivot_row = 0
    for j in range(cols):
        if pivot_row < rows:
            i = pivot_row
            while i < rows and mat[i, j] == 0:
                i += 1
            
            if i < rows:
                mat[[pivot_row, i]] = mat[[i, pivot_row]]
                for k in range(rows):
                    if k != pivot_row and mat[k, j] == 1:
                        mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                pivot_row += 1
    return pivot_row

def solve_expected_rounds():
    """
    Calculates the expected number of rounds for the system to return to its initial state.
    """
    # Define influence sets (using 0-based indexing for matrix)
    influence_sets = {
        0: {1, 3, 5, 6},   # Person 1
        1: {2, 4, 5, 7},   # Person 2
        2: {3, 5},         # Person 3
        3: {4},            # Person 4
        4: {5, 7},         # Person 5
        5: {6},            # Person 6
        6: {7},            # Person 7
        7: {}              # Person 8
    }

    # Construct the influence matrix M
    # M[j, i] = 1 if person i+1 influences person j+1
    M = np.zeros((8, 8), dtype=int)
    for influencer, influenced_set in influence_sets.items():
        for influenced in influenced_set:
            M[influenced, influencer] = 1

    # Compute powers of M
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2

    # Calculate ranks of M, M^2, M^4
    rank_M1 = rank_F2(M)
    rank_M2 = rank_F2(M2)
    rank_M4 = rank_F2(M4)

    # Calculate the size of the kernel for each matrix
    # |ker(M^k)| = 2^(8 - rank(M^k))
    size_ker_M1 = 2**(8 - rank_M1)
    size_ker_M2 = 2**(8 - rank_M2)
    size_ker_M4 = 2**(8 - rank_M4)
    size_ker_M8 = 2**8 # Since M^5 = 0, M^8 = 0, so rank(M^8) = 0

    # Calculate the number of states for each order
    N1 = size_ker_M1
    N2 = size_ker_M2 - size_ker_M1
    N4 = size_ker_M4 - size_ker_M2
    N8 = size_ker_M8 - size_ker_M4
    
    orders = [1, 2, 4, 8]
    counts = [N1, N2, N4, N8]
    
    print("Number of states for each possible order:")
    for order, count in zip(orders, counts):
        print(f"Order {order}: {count} states")
    print("-" * 30)

    # Calculate the total sum for the expectation formula
    total_sum = (1 * N1) + (2 * N2) + (4 * N4) + (8 * N8)
    
    # Calculate the expected value
    expected_R = total_sum / 256.0
    
    print("Calculation of the expected value E[R]:")
    print(f"E[R] = (1/256) * (1 * {N1} + 2 * {N2} + 4 * {N4} + 8 * {N8})")
    print(f"E[R] = (1/256) * ({1 * N1} + {2 * N2} + {4 * N4} + {8 * N8})")
    print(f"E[R] = (1/256) * {total_sum}")
    print(f"E[R] = {total_sum / 256}")
    print(f"E[R] rounded to 2 decimal places: {expected_R:.2f}")

    return expected_R

if __name__ == "__main__":
    solve_expected_rounds()
    print("\n<<<6.86>>>")
