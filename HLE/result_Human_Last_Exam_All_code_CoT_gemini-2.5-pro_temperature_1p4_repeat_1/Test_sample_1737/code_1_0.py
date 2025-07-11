import numpy as np
from itertools import combinations_with_replacement

def solve_for_su_n(N):
    """
    Calculates the number of distinct non-zero values for the symmetric
    structure constants d_ijk of SU(N).
    
    Args:
        N (int): The dimension of the special unitary group SU(N).

    Returns:
        int: The number of distinct non-zero d_ijk values.
    """
    if not isinstance(N, int) or N < 2:
        print("Please provide an integer N >= 2.")
        return 0

    # List to store the N^2 - 1 generator matrices
    gens = []

    # 1. Symmetric off-diagonal generators
    # There are N*(N-1)/2 of these.
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = 0.5
            mat[k, j] = 0.5
            # Normalization Tr(T^2) = 1/2 is satisfied by this construction.
            gens.append(mat)

    # 2. Antisymmetric off-diagonal generators
    # There are N*(N-1)/2 of these.
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = -0.5j
            mat[k, j] = 0.5j
            # Normalization Tr(T^2) = 1/2 is satisfied by this construction.
            gens.append(mat)

    # 3. Diagonal generators
    # There are N-1 of these.
    for l in range(1, N):
        mat = np.zeros((N, N), dtype=complex)
        # The normalization factor
        norm = np.sqrt(2 * l * (l + 1))
        for i in range(l):
            mat[i, i] = 1.0
        mat[l, l] = -l
        mat /= norm
        gens.append(mat)

    dim = N**2 - 1
    
    # Set to store unique, non-zero d_ijk values
    unique_d_values = set()
    
    indices = range(dim)
    
    # Iterate over all combinations of three generators (i, j, k) with i <= j <= k
    for i, j, k in combinations_with_replacement(indices, 3):
        Ti, Tj, Tk = gens[i], gens[j], gens[k]
        
        # Calculate d_ijk = 2 * Tr(T_i{T_j, T_k}) = 2 * Tr(T_i T_j T_k + T_j T_i T_k)
        # The trace should be real, but we use np.real to discard negligible imaginary parts
        # due to floating point inaccuracies.
        trace_val = np.trace(Ti @ Tj @ Tk + Tj @ Ti @ Tk)
        d_val = 2 * np.real(trace_val)
        
        # Round the result to handle floating point precision issues
        d_rounded = round(d_val, 9)
        
        # Add the value to our set if it's non-zero
        if abs(d_rounded) > 1e-10:
            unique_d_values.add(d_rounded)
            
    return len(unique_d_values)

if __name__ == '__main__':
    # Set the value of N for the SU(N) group.
    # For SU(2), the result is 0.
    # For SU(3), the result is 5.
    # For SU(4), the result is 9.
    N = 3
    
    num_values = solve_for_su_n(N)
    print(f"For SU({N}), the number of different numerical values for non-zero d_ijk is: {num_values}")
