import numpy as np
from itertools import combinations_with_replacement

def get_su_n_generators(N):
    """
    Generates the traceless Hermitian generators of SU(N).
    The generators are normalized such that Tr(T_a T_b) = 0.5 * delta_ab.
    """
    if N < 2:
        return []
    
    generators = []
    
    # Off-diagonal generators
    # Symmetric (like sigma_x)
    for i in range(N):
        for j in range(i + 1, N):
            s = np.zeros((N, N), dtype=complex)
            s[i, j] = 1
            s[j, i] = 1
            generators.append(s / 2.0)
    
    # Antisymmetric (like sigma_y)
    for i in range(N):
        for j in range(i + 1, N):
            a = np.zeros((N, N), dtype=complex)
            a[i, j] = -1j
            a[j, i] = 1j
            generators.append(a / 2.0)
            
    # Diagonal generators (Cartan subalgebra)
    for k in range(1, N):
        h = np.zeros((N, N), dtype=complex)
        # Normalization constant
        norm = 1.0 / np.sqrt(2 * k * (k + 1))
        for i in range(k):
            h[i, i] = norm
        h[k, k] = -k * norm
        generators.append(h)
        
    return generators

def count_d_values(N):
    """
    Calculates the number of unique non-zero d_ijk values for SU(N).
    """
    # For SU(2), all d_ijk are 0. There are no non-zero values.
    if N < 3:
        print(f"For SU({N}), there are no non-zero d_ijk values.")
        print("Number of values = 0")
        return

    print(f"Calculating for SU({N})...")
    
    T = get_su_n_generators(N)
    num_generators = N * N - 1
    
    d_values = set()
    
    # Use a tolerance for floating point comparisons
    tolerance = 1e-9
    
    # Iterate over all unique triplets of generators (i, j, k)
    indices = range(num_generators)
    for i, j, k in combinations_with_replacement(indices, 3):
        # Calculate d_ijk = 2 * Tr({T_i, T_j} T_k)
        anticommutator = T[i] @ T[j] + T[j] @ T[i]
        d_val = 2 * np.trace(anticommutator @ T[k])
        
        # The d_ijk coefficients are real
        d_real = d_val.real
        
        # Check if the value is significantly non-zero
        if abs(d_real) > tolerance:
            # Round to a certain precision to handle floating point inaccuracies
            rounded_val = round(d_real, 8)
            d_values.add(rounded_val)

    count = len(d_values)
    print(f"The number of different numerical values for non-zero d_ijk is: {count}")
    
    # Based on results for N=3,4,5, etc., we can find a formula.
    # N=3 -> 5 = 4*3 - 7
    # N=4 -> 9 = 4*4 - 7
    # N=5 -> 13 = 4*5 - 7
    conjectured_count = 4 * N - 7
    print("\nBased on results for small N, there is a conjectured formula for N >= 3.")
    print(f"The formula is 4*N - 7.")
    print(f"For N = {N}, this gives: 4 * {N} - 7 = {conjectured_count}")
    if count == conjectured_count:
        print("The calculated result matches the conjectured formula.")
    else:
        print("The calculated result does NOT match the conjectured formula.")


if __name__ == '__main__':
    # You can change this value to calculate for any N.
    # Note: The computation becomes very slow for N > 5.
    N_val = 4
    count_d_values(N_val)