import numpy as np

def solve_su_n_structure_constants():
    """
    Calculates and counts the number of different numerical values for the 
    non-zero totally symmetric structure constants d_ijk of the SU(N) group.
    """
    # Set the value of N for the SU(N) group.
    # You can change this value to analyze other SU(N) groups (for N >= 2).
    N = 4

    print(f"Analyzing the symmetric structure constants for SU({N}).")

    # Handle the special case of SU(2), where all d_ijk are known to be zero.
    if N == 2:
        print("\nFor SU(2), the anti-commutator of any two generators is a multiple of the identity matrix.")
        print("{T_i, T_j} = (1/2) * delta_ij * I")
        print("This implies that all structure constants d_ijk are zero.")
        print("\nNumber of different non-zero values: 0")
        return 0

    # Step 1: Generate a basis of SU(N) generators T_a
    # The generators are normalized such that Tr(T_a T_b) = (1/2) delta_ab.
    generators = []
    
    # Part 1: Off-diagonal generators (N*(N-1) of them)
    for j in range(N):
        for k in range(j + 1, N):
            # Symmetric generator
            s_jk = np.zeros((N, N), dtype=complex)
            s_jk[j, k] = 1
            s_jk[k, j] = 1
            generators.append(s_jk / 2)

            # Anti-symmetric generator
            a_jk = np.zeros((N, N), dtype=complex)
            a_jk[j, k] = -1j
            a_jk[k, j] = 1j
            generators.append(a_jk / 2)

    # Part 2: Diagonal generators (N-1 of them)
    for m in range(1, N):
        h_m = np.zeros((N, N), dtype=complex)
        # This prefactor ensures Tr(H_m^2) = 2, so Tr(T^2) will be 1/2.
        prefactor = np.sqrt(2 / (m * (m + 1)))
        for i in range(m):
            h_m[i, i] = 1
        h_m[m, m] = -m
        h_m *= prefactor
        generators.append(h_m / 2)

    num_generators = len(generators)
    print(f"\nConstructed a basis of {num_generators} generators for SU({N}).")

    # Step 2: Calculate all unique non-zero d_ijk values
    d_values = set()
    # Tolerance for floating point comparisons
    tolerance = 1e-9
    rounding_decimals = 8

    print("Calculating structure constants d_ijk = 2 * Tr({T_i, T_j} T_k)...")
    # Iterate over i, j, k with i <= j <= k to avoid permutations
    for i in range(num_generators):
        for j in range(i, num_generators):
            for k in range(j, num_generators):
                T_i, T_j, T_k = generators[i], generators[j], generators[k]
                
                # The d_ijk are real, so we take the real part of the trace result
                # to discard any small imaginary part arising from numerical noise.
                anticommutator = T_i @ T_j + T_j @ T_i
                d_val = 2 * np.real(np.trace(anticommutator @ T_k))
                
                # Check if the value is non-zero within the given tolerance
                if abs(d_val) > tolerance:
                    # Round the value to handle minor floating point variations
                    rounded_val = round(d_val, rounding_decimals)
                    # Add to set to automatically handle uniqueness
                    if abs(rounded_val) > tolerance: # Check again after rounding
                        d_values.add(rounded_val)
    
    # Step 3: Print the results
    num_distinct_values = len(d_values)
    sorted_values = sorted(list(d_values))

    print("\n--- Results ---")
    print(f"The number of different numerical values for the non-zero d_ijk is: {num_distinct_values}")
    
    if num_distinct_values > 0:
        print("The distinct values are:")
        for val in sorted_values:
            print(f"  {val:.8f}")
            
    return num_distinct_values

if __name__ == '__main__':
    solve_su_n_structure_constants()
