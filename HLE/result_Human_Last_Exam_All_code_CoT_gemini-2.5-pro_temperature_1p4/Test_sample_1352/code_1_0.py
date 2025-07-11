import numpy as np

def rank_F2(matrix):
    """
    Calculates the rank of a binary matrix over the field F_2.
    """
    mat = matrix.copy().astype(int)
    rows, cols = mat.shape
    rank = 0
    pivot_row = 0
    for j in range(cols):
        if pivot_row < rows:
            i = pivot_row
            while i < rows and mat[i, j] == 0:
                i += 1
            if i < rows:
                # Swap rows to bring pivot to the front
                mat[[pivot_row, i]] = mat[[i, pivot_row]]
                # Eliminate other 1s in the same column
                for k in range(rows):
                    if k != pivot_row and mat[k, j] == 1:
                        mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                pivot_row += 1
    return pivot_row

def solve():
    """
    Solves the problem by calculating the expected number of rounds.
    """
    # Influence sets (1-based index)
    influence_sets = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6},
        4: {5}, 5: {6, 8}, 6: {7},
        7: {8}, 8: {}
    }

    # Build the matrix A (0-based index for numpy)
    # A[j, i] = 1 if person i+1 influences person j+1
    A = np.zeros((8, 8), dtype=int)
    for i in range(1, 9):
        for j in influence_sets[i]:
            A[j - 1, i - 1] = 1

    # Compute powers of A modulo 2
    A2 = (A @ A) % 2
    A4 = (A2 @ A2) % 2
    A8 = (A4 @ A4) % 2 # This will be the zero matrix

    # Calculate ranks over F_2
    rank_A = rank_F2(A)
    rank_A2 = rank_F2(A2)
    rank_A4 = rank_F2(A4)
    # rank_A8 is 0 as A^8 is the zero matrix

    # Calculate the size of the kernels
    # |Ker(B)| = 2^(dim - rank(B))
    num_states = 2**8
    ker_A_size = 2**(8 - rank_A)
    ker_A2_size = 2**(8 - rank_A2)
    ker_A4_size = 2**(8 - rank_A4)
    ker_A8_size = num_states

    # n_d is the number of vectors with order d
    # Sum of n_k for k|d is |Ker(A^d)|
    n1 = ker_A_size
    n2 = ker_A2_size - n1
    n4 = ker_A4_size - (n1 + n2)
    n8 = ker_A8_size - (n1 + n2 + n4)

    # Calculate the expected value E[R]
    # E[R] = sum(d * n_d) / total_states
    sum_of_orders = (1 * n1) + (2 * n2) + (4 * n4) + (8 * n8)
    expected_R = sum_of_orders / num_states
    
    # Print the result step-by-step
    print("The number of states with each order are:")
    print(f"Order 1: n1 = {n1}")
    print(f"Order 2: n2 = {n2}")
    print(f"Order 4: n4 = {n4}")
    print(f"Order 8: n8 = {n8}")
    print("\nThe expected value E[R] is the average of all orders:")
    
    # Using variables to construct the formula string
    formula = f"(1*{n1} + 2*{n2} + 4*{n4} + 8*{n8}) / {num_states}"
    calc1 = f"({1*n1} + {2*n2} + {4*n4} + {8*n8}) / {num_states}"
    calc2 = f"{sum_of_orders} / {num_states}"
    
    print(f"E[R] = {formula}")
    print(f"E[R] = {calc1}")
    print(f"E[R] = {calc2}")
    print(f"E[R] = {expected_R}")
    print(f"E[R] rounded to 2 decimal places = {expected_R:.2f}")

solve()
<<<6.93>>>