import numpy as np
from itertools import combinations_with_replacement

def get_su_n_generators(N):
    """
    Generates a basis of traceless Hermitian generators for SU(N),
    normalized such that Tr(T_a T_b) = 0.5 * delta_ab.
    This uses the generalized Gell-Mann matrix construction.
    """
    if N < 2:
        return []

    num_gens = N * N - 1
    generators = []
    
    # 1. Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            gen[j, k] = 1.0
            gen[k, j] = 1.0
            # Normalize to Tr(T^2) = 1/2. Tr(gen^2) is 2. So scale by 1/2.
            generators.append(gen / 2.0)

    # 2. Antisymmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            gen[j, k] = -1.0j
            gen[k, j] = 1.0j
            # Normalize to Tr(T^2) = 1/2. Tr(gen^2) is 2. So scale by 1/2.
            generators.append(gen / 2.0)

    # 3. Diagonal generators
    for l in range(1, N):
        gen = np.zeros((N, N), dtype=complex)
        for i in range(l):
            gen[i, i] = 1.0
        gen[l, l] = -l
        # Normalize to Tr(T^2) = 1/2. Tr(gen^2) is l^2 + l. So scale by 1/sqrt(2*(l^2+l))
        norm = np.sqrt(2 * l * (l + 1))
        generators.append(gen / norm)
        
    # Safety check
    if len(generators) != num_gens:
        raise ValueError(f"Generated {len(generators)} generators, but expected {num_gens}")

    return generators

def solve_for_n(N):
    """
    Calculates the number of distinct non-zero values for d_ijk of SU(N).
    """
    print(f"Calculating for SU({N})...")
    # For SU(2), all d_ijk are zero. The number of non-zero values is 0.
    if N < 3:
        print("For SU(N) with N < 3, all d_ijk coefficients are zero.")
        print("\nNumber of distinct non-zero values: 0")
        return 0

    try:
        generators = get_su_n_generators(N)
    except ValueError as e:
        print(f"Error: {e}")
        return
        
    num_gens = len(generators)
    
    # Tolerance for floating point comparisons
    TOL = 1e-9
    distinct_values = []

    indices = range(num_gens)
    # Iterate through all combinations i <= j <= k because d_ijk is totally symmetric
    for i, j, k in combinations_with_replacement(indices, 3):
        Ti = generators[i]
        Tj = generators[j]
        Tk = generators[k]

        # Calculate d_ijk = 2 * Tr({T_i, T_j} T_k)
        anticomm = Ti @ Tj + Tj @ Ti
        d_val = 2 * np.trace(anticomm @ Tk)

        # d_ijk must be real since generators are Hermitian. We take the real part.
        val = d_val.real
        
        # If the value is non-zero (outside the tolerance range)
        if abs(val) > TOL:
            # Check if this value (or a very close one) is already in our list
            is_new = True
            for existing_val in distinct_values:
                if abs(val - existing_val) < TOL:
                    is_new = False
                    break
            
            if is_new:
                distinct_values.append(val)
    
    print("\nThe distinct non-zero numerical values for d_ijk are:")
    distinct_values.sort()
    for v in distinct_values:
        print(f"{v:.8f}")
                
    count = len(distinct_values)
    print(f"\nTotal number of distinct non-zero values: {count}")
    return count

if __name__ == '__main__':
    # Set the value of N for the SU(N) group you are interested in.
    # For example, N=3 for the group of the strong interaction, N=4, etc.
    N = 4
    
    solve_for_n(N)