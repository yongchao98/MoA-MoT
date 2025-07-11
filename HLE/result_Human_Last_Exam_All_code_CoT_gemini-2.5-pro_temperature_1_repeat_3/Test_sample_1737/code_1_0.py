import numpy as np
import math
from itertools import combinations_with_replacement

def count_distinct_d_values(N):
    """
    Calculates the number of distinct non-zero values for the symmetric
    structure constants d_ijk of SU(N).

    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        int: The number of distinct non-zero d_ijk values.
    """
    if not isinstance(N, int) or N < 2:
        print("Error: N must be an integer greater than or equal to 2.")
        return 0

    # Step 1: Generate the SU(N) generators (generalized Gell-Mann matrices)
    dim = N * N - 1
    generators = []

    # Type 1: Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=np.complex128)
            mat[j, k] = 1.0
            mat[k, j] = 1.0
            # Normalization: Tr(T^2) = 1/2. For this type, Tr(mat^2) = 2.
            # So we normalize by 1/2.
            generators.append(mat / 2.0)

    # Type 2: Antisymmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=np.complex128)
            mat[j, k] = -1.0j
            mat[k, j] = 1.0j
            # Normalization: Tr(T^2) = 1/2. For this type, Tr(mat^2) = 2.
            # So we normalize by 1/2.
            generators.append(mat / 2.0)

    # Type 3: Diagonal generators
    for l in range(1, N):
        mat = np.zeros((N, N), dtype=np.complex128)
        for i in range(l):
            mat[i, i] = 1.0
        mat[l, l] = -l
        # Normalization: Tr(T^2) = 1/2. For this type, Tr(mat^2) = l^2 + l.
        # So we normalize by 1/sqrt(2 * (l^2+l)).
        norm = math.sqrt(2 * l * (l + 1))
        generators.append(mat / norm)

    # Step 2 & 3: Calculate d_ijk and count unique non-zero values
    d_values = set()
    # To handle floating point inaccuracies, we round the results
    precision = 12

    # Iterate over all combinations of indices i, j, k with replacement
    indices = range(dim)
    for i, j, k in combinations_with_replacement(indices, 3):
        # The formula is d_ijk = 4 * Re(Tr(T_i T_j T_k))
        # np.dot is matrix multiplication for 2D arrays
        trace_val = np.trace(generators[i] @ generators[j] @ generators[k])
        d_val = 4 * np.real(trace_val)

        # Add the value to the set if it's not close to zero
        if not np.isclose(d_val, 0):
            d_values.add(round(d_val, precision))

    # Step 4: Output the result
    num_distinct_values = len(d_values)
    print(f"For SU({N}), the number of different numerical values for non-zero d_ijk is = {num_distinct_values}")
    
    return num_distinct_values

if __name__ == '__main__':
    # You can change this value to calculate the result for a different SU(N)
    N = 4
    
    count_distinct_d_values(N)