import numpy as np
from itertools import combinations_with_replacement

def get_sun_generators(N):
    """
    Generates the generalized Gell-Mann matrices for SU(N).
    These are a set of N^2-1 traceless Hermitian matrices lambda_a
    normalized such that Tr(lambda_a lambda_b) = 2 delta_ab.
    """
    if N < 2:
        return []

    num_gens = N**2 - 1
    generators = []

    # Type 1: Symmetric off-diagonal matrices
    # N(N-1)/2 of these
    for j in range(N):
        for k in range(j + 1, N):
            matrix = np.zeros((N, N), dtype=complex)
            matrix[j, k] = 1
            matrix[k, j] = 1
            generators.append(matrix)

    # Type 2: Antisymmetric off-diagonal matrices
    # N(N-1)/2 of these
    for j in range(N):
        for k in range(j + 1, N):
            matrix = np.zeros((N, N), dtype=complex)
            matrix[j, k] = -1j
            matrix[k, j] = 1j
            generators.append(matrix)

    # Type 3: Diagonal matrices
    # N-1 of these
    for l in range(1, N):
        matrix = np.zeros((N, N), dtype=complex)
        # Normalization factor
        norm = np.sqrt(2 / (l * (l + 1)))
        for i in range(l):
            matrix[i, i] = norm
        matrix[l, l] = -l * norm
        generators.append(matrix)
        
    if len(generators) != num_gens:
        raise ValueError(f"Generated {len(generators)} generators, but expected {num_gens} for SU({N}).")

    return generators

def calculate_d_ijk_values(N):
    """
    Calculates the set of unique non-zero d_ijk values for SU(N).
    """
    if N < 2:
        print(f"SU({N}) is not well-defined in this context. There are no generators.")
        return set()
    
    # For SU(2), the anti-commutator is proportional to the identity matrix,
    # which means all d_ijk coefficients are zero.
    if N == 2:
        return set()

    generators = get_sun_generators(N)
    num_gens = len(generators)
    
    d_values = set()
    
    # Iterate over combinations with replacement for (i, j, k) to avoid redundant calculations.
    indices = range(num_gens)
    for i, j, k in combinations_with_replacement(indices, 3):
        l_i = generators[i]
        l_j = generators[j]
        l_k = generators[k]
        
        # d_ijk = (1/4) * Tr({lambda_i, lambda_j} * lambda_k)
        anticommutator = l_i @ l_j + l_j @ l_i
        d_ijk_val = 0.25 * np.trace(anticommutator @ l_k)
        
        d_val_real = np.real(d_ijk_val)

        if not np.isclose(d_val_real, 0):
            # Round values to handle floating point inaccuracies before adding to the set.
            rounded_val = round(d_val_real, 8)
            d_values.add(rounded_val)

    return d_values

def main():
    # Set the desired value of N for the SU(N) group.
    # For example, N=3 for the group of Quantum Chromodynamics.
    # You can change this value to explore other groups like SU(4), SU(5), etc.
    N = 3

    print(f"Analyzing the symmetric structure constants d_ijk for SU({N}).")
    
    unique_d_values = calculate_d_ijk_values(N)
    
    num_unique_values = len(unique_d_values)
    
    print("\nThe unique non-zero numerical values for d_ijk are:")
    if num_unique_values > 0:
        # Print sorted values for better readability
        print(sorted(list(unique_d_values)))
    else:
        print("None (all d_ijk are zero).")
        
    print(f"\nFor SU({N}), the number of different non-zero numerical values for d_ijk is: {num_unique_values}")


if __name__ == '__main__':
    main()