import numpy as np
from itertools import combinations_with_replacement

def count_unique_d_values_for_sun(N):
    """
    Calculates the number of unique non-zero values for the symmetric structure
    constants d_ijk of SU(N).

    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        None. Prints the unique values and their count.
    """
    if not isinstance(N, int) or N < 2:
        print("Please provide an integer N >= 2.")
        return

    # --- Step 1: Generate the SU(N) generators (generalized Gell-Mann matrices) ---
    num_generators = N**2 - 1
    generators = []

    # Helper to create an elementary matrix E_jk
    def E(j, k, size):
        mat = np.zeros((size, size), dtype=complex)
        mat[j, k] = 1
        return mat

    # Off-diagonal generators (indices are 0-based)
    for j in range(N):
        for k in range(j + 1, N):
            # Symmetric generator
            T_s = 0.5 * (E(j, k, N) + E(k, j, N))
            generators.append(T_s)
            
            # Anti-symmetric generator
            T_a = -0.5j * (E(j, k, N) - E(k, j, N))
            generators.append(T_a)

    # Diagonal generators
    for l_idx in range(1, N): # l runs from 1 to N-1
        norm_factor = 1.0 / np.sqrt(2 * l_idx * (l_idx + 1))
        T_d = np.zeros((N, N), dtype=complex)
        for j in range(l_idx):
            T_d[j, j] = 1
        T_d[l_idx, l_idx] = -l_idx
        T_d *= norm_factor
        generators.append(T_d)
    
    # --- Step 2: Calculate d_ijk values and find unique ones ---
    unique_d_values = set()
    # Floating point precision for checking zero and for rounding
    TOLERANCE = 1e-9
    PRECISION = 8

    # Iterate over all combinations of three generators with replacement
    indices = range(num_generators)
    for i, j, k in combinations_with_replacement(indices, 3):
        T_i, T_j, T_k = generators[i], generators[j], generators[k]
        
        # d_ijk = 4 * Tr(T_i * T_j * T_k)
        d_val = 4 * np.trace(T_i @ T_j @ T_k)
        
        # The result must be real. Take real part to discard numerical noise.
        d_real = np.real(d_val)
        
        if abs(d_real) > TOLERANCE:
            d_rounded = round(d_real, PRECISION)
            # Add to set to ensure uniqueness
            unique_d_values.add(d_rounded)

    # --- Step 3: Print the results ---
    print(f"For SU({N}):")
    if not unique_d_values:
        print("All d_ijk coefficients are zero.")
        print("The number of different non-zero numerical values is: 0")
    else:
        # Sort for clearer presentation
        sorted_values = sorted(list(unique_d_values))
        print("The unique non-zero values for d_ijk are:")
        print(sorted_values)
        print(f"The number of different non-zero numerical values is: {len(unique_d_values)}")

if __name__ == '__main__':
    try:
        N_input = int(input("Enter the value of N for SU(N): "))
        count_unique_d_values_for_sun(N_input)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
