import numpy as np
from itertools import combinations_with_replacement

def get_su_n_generators(N):
    """
    Generates a basis of traceless Hermitian generators for SU(N)
    normalized such that Tr(T_a T_b) = (1/2) delta_ab.
    This basis is a generalization of the Gell-Mann matrices.
    """
    if N < 2:
        return []

    generators = []
    
    # Off-diagonal generators
    # There are N*(N-1)/2 symmetric and N*(N-1)/2 anti-symmetric ones.
    for i in range(N):
        for j in range(i + 1, N):
            # Symmetric generator
            s_ij = np.zeros((N, N), dtype=complex)
            s_ij[i, j] = 1
            s_ij[j, i] = 1
            generators.append(s_ij / 2)
            
            # Anti-symmetric generator
            a_ij = np.zeros((N, N), dtype=complex)
            a_ij[i, j] = -1j
            a_ij[j, i] = 1j
            generators.append(a_ij / 2)
            
    # Diagonal generators
    # There are N-1 of them.
    for m in range(1, N):
        h_m = np.zeros((N, N), dtype=complex)
        for k in range(m):
            h_m[k, k] = 1
        h_m[m, m] = -m
        norm = np.sqrt(2 * m * (m + 1))
        generators.append(h_m / norm)
        
    return generators

def find_d_ijk_values(N):
    """
    Calculates the set of unique non-zero d_ijk values for SU(N).
    """
    if N < 2:
        print("SU(N) is trivial or has no d_ijk constants for N < 2.")
        return 0, set()
    
    # For SU(2), all d_ijk are zero.
    if N == 2:
        return 0, set()

    generators = get_su_n_generators(N)
    num_gens = len(generators)
    
    d_values = set()
    
    # Iterate over all combinations of indices i, j, k with replacement
    indices = range(num_gens)
    for i, j, k in combinations_with_replacement(indices, 3):
        Ti = generators[i]
        Tj = generators[j]
        Tk = generators[k]
        
        # d_ijk = 2 * Tr({Ti, Tj} Tk) where {Ti, Tj} = Ti @ Tj + Tj @ Ti
        anticommutator = Ti @ Tj + Tj @ Ti
        d = 2 * np.trace(anticommutator @ Tk)

        # The structure constants d_ijk are real. Imaginary part should be zero.
        d_real = d.real
        
        # Round to handle floating point inaccuracies
        d_rounded = round(d_real, 10)
        
        # Add to the set if it's non-zero
        if abs(d_rounded) > 1e-9:
            d_values.add(d_rounded)
            
    return len(d_values), sorted(list(d_values))

def main():
    """
    Main function to get user input and run the calculation.
    """
    try:
        n_input = input("Please enter the value of N for SU(N): ")
        N = int(n_input)
        if N < 1:
            raise ValueError("N must be a positive integer.")
    except (ValueError, TypeError):
        print(f"Error: Invalid input '{n_input}'. Please enter an integer greater than 0.")
        return

    print(f"\nCalculating the number of distinct non-zero d_ijk values for SU({N})...")
    
    count, values = find_d_ijk_values(N)
    
    print(f"\nFor SU({N}), there are {count} different numerical values for the non-zero d_ijk.")
    if count > 0:
        print("The values are:")
        print(values)

if __name__ == "__main__":
    main()