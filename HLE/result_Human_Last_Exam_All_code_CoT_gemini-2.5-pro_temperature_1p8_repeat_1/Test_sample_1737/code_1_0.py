import numpy as np
import itertools

def count_distinct_d_values(N: int):
    """
    Calculates the number of distinct numerical values for the non-zero
    symmetric structure constants d_ijk of SU(N).

    Args:
        N: The dimension of the special unitary group SU(N).
    """
    if not isinstance(N, int) or N < 2:
        print("Error: N must be an integer greater than or equal to 2.")
        return

    # For SU(2), all d_ijk are 0.
    if N == 2:
        print(0)
        return

    # Precision for rounding to handle floating point comparisons
    precision = 10
    
    # Use a set to store unique values
    # The off-diagonal-only types give values of +/- 0.5 for N >= 3
    distinct_values = {0.5, -0.5}

    # Generate the N-1 diagonal lambda matrices
    diag_gens = []
    for m in range(1, N):
        c = np.sqrt(2.0 / (m * (m + 1)))
        diag = np.zeros(N)
        diag[:m] = 1.0
        diag[m] = -m
        diag_gens.append(c * np.diag(diag))

    # Calculate values from d(D, S, S) and d(D, A, A) types
    # The formula is the same for S(symmetric) and A(anti-symmetric) types
    for m in range(N - 1):
        gen_m = diag_gens[m]
        for j in range(N):
            for k in range(j + 1, N):
                # Python indices are 0-based
                val = 0.5 * (gen_m[j, j] + gen_m[k, k])
                if not np.isclose(val, 0):
                    distinct_values.add(round(val, precision))

    # Calculate values from d(D, D, D) types
    indices = range(N - 1)
    # Consider combinations with replacement of 3 indices for the diagonal generators
    for m1, m2, m3 in itertools.combinations_with_replacement(indices, 3):
        gen1 = diag_gens[m1]
        gen2 = diag_gens[m2]
        gen3 = diag_gens[m3]
        
        # d_ijk = 0.5 * Tr(L_i L_j L_k) for diagonal commuting generators
        val = 0.5 * np.trace(gen1 @ gen2 @ gen3)
        if not np.isclose(val, 0):
            distinct_values.add(round(val, precision))
    
    print(len(distinct_values))


if __name__ == '__main__':
    # You can change the value of N here to get the answer for different groups.
    # For example, N=3, N=4, etc.
    try:
        N_val = int(input("Enter the value of N for SU(N) (must be >= 2): "))
        count_distinct_d_values(N_val)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
