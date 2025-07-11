import numpy as np

def rank_F2(A):
    """
    Computes the rank of a matrix over the finite field F_2.
    """
    m, n = A.shape
    A_copy = A.copy().astype(int)
    rank = 0
    pivot_row = 0
    for j in range(n):  # Iterate over columns
        if pivot_row < m:
            pivot = pivot_row
            while pivot < m and A_copy[pivot, j] == 0:
                pivot += 1
            
            if pivot < m:
                # Swap rows to bring pivot to the top of the submatrix
                A_copy[[pivot_row, pivot]] = A_copy[[pivot, pivot_row]]
                
                # Eliminate other 1s in the same column
                for i in range(m):
                    if i != pivot_row and A_copy[i, j] == 1:
                        A_copy[i, :] = (A_copy[i, :] + A_copy[pivot_row, :]) % 2
                pivot_row += 1
    return pivot_row

def solve():
    """
    Solves the problem by calculating the expected number of rounds.
    """
    # Define influence sets (1-indexed)
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: {}
    }

    # Create the influence matrix M (0-indexed)
    # M_ji = 1 if person i influences person j
    M = np.zeros((8, 8), dtype=int)
    for i in range(1, 9):
        for j in influence_sets[i]:
            M[j - 1, i - 1] = 1

    # The state transition is s_new = (I + M) * s
    # The order of a vector s is the smallest r >= 1 such that (I+M)^r * s = s
    # This is equivalent to ((I+M)^r - I) * s = 0
    # Over F_2, (I+M)^(2^k) - I = M^(2^k)
    # We need to find the number of vectors in the null spaces of M, M^2, M^4, etc.

    # Compute powers of M
    M2 = (M @ M) % 2
    M4 = (M2 @ M2) % 2
    M8 = (M4 @ M4) % 2 # This will be the zero matrix

    # Compute ranks
    rank_M = rank_F2(M)
    rank_M2 = rank_F2(M2)
    rank_M4 = rank_F2(M4)
    rank_M8 = rank_F2(M8)

    # The number of vectors v with order dividing d=2^k is |Null(M^d)| = 2^(8 - rank(M^d))
    num_vectors_order_div_1 = 2**(8 - rank_M)
    num_vectors_order_div_2 = 2**(8 - rank_M2)
    num_vectors_order_div_4 = 2**(8 - rank_M4)
    num_vectors_order_div_8 = 2**(8 - rank_M8)

    # Number of vectors with exact order d
    N1 = num_vectors_order_div_1
    N2 = num_vectors_order_div_2 - num_vectors_order_div_1
    N4 = num_vectors_order_div_4 - num_vectors_order_div_2
    N8 = num_vectors_order_div_8 - num_vectors_order_div_4

    # Calculate the expected value E[R]
    total_sum_of_orders = (1 * N1 + 2 * N2 + 4 * N4 + 8 * N8)
    expected_R = total_sum_of_orders / 256

    # Print the results
    print("The number of initial states for each possible return time R:")
    print(f"R = 1: {N1} states")
    print(f"R = 2: {N2} states")
    print(f"R = 4: {N4} states")
    print(f"R = 8: {N8} states")
    print("\nThe calculation for the expected value E[R] is:")
    print(f"E[R] = (1 * {N1} + 2 * {N2} + 4 * {N4} + 8 * {N8}) / 256")
    print(f"E[R] = ({1 * N1} + {2 * N2} + {4 * N4} + {8 * N8}) / 256")
    print(f"E[R] = {total_sum_of_orders} / 256")
    print(f"E[R] = {expected_R}")
    print(f"\nRounded to 2 decimal places, the expected value is: {expected_R:.2f}")

solve()
<<<7.46>>>