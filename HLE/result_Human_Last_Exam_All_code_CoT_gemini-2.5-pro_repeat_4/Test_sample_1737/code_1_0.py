import numpy as np

def count_distinct_d_values_for_sun(N):
    """
    Calculates the number of distinct non-zero values for the symmetric 
    structure constants d_ijk of SU(N).

    The totally symmetric structure constants d_ijk are defined by the
    anti-commutation relation of the group generators T_a:
    {T_i, T_j} = (1/N) * delta_ij * I + d_ijk * T_k

    From this, one can derive:
    d_ijk = 2 * Tr({T_i, T_j} T_k)
    
    This script implements this calculation.
    """
    print(f"Analyzing for SU(N) with N = {N}")

    if not isinstance(N, int) or N < 2:
        print("Error: N must be an integer greater than or equal to 2.")
        return

    # For SU(2), the d_ijk coefficients are all zero.
    if N == 2:
        print("For SU(2), all d_ijk coefficients are zero.")
        print("Number of distinct non-zero values: 0")
        print("The distinct values are: []")
        return

    # Step 1: Generate the SU(N) generators (generalized Gell-Mann matrices)
    # They are normalized such that Tr(T_a T_b) = (1/2) * delta_ab
    generators = []
    
    # Off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            # Symmetric
            s_matrix = np.zeros((N, N), dtype=complex)
            s_matrix[j, k] = 1
            s_matrix[k, j] = 1
            generators.append(s_matrix / 2.0)
            
            # Anti-symmetric
            a_matrix = np.zeros((N, N), dtype=complex)
            a_matrix[j, k] = -1j
            a_matrix[k, j] = 1j
            generators.append(a_matrix / 2.0)

    # Diagonal generators
    for l in range(1, N):
        d_matrix = np.zeros((N, N), dtype=complex)
        diag_elements = np.ones(l)
        diag_elements = np.append(diag_elements, -l)
        diag_elements = np.append(diag_elements, np.zeros(N - l - 1))
        norm_factor = 1.0 / np.sqrt(2 * l * (l + 1))
        np.fill_diagonal(d_matrix, diag_elements * norm_factor)
        generators.append(d_matrix)

    num_gen = len(generators)
    d_values = set()
    
    # Tolerance for floating point comparisons
    TOL = 1e-9

    # Step 2: Iterate and calculate d_ijk values
    # d_ijk is totally symmetric, so we only iterate over i <= j <= k
    for i in range(num_gen):
        for j in range(i, num_gen):
            # Precompute the anti-commutator
            anti_comm = generators[i] @ generators[j] + generators[j] @ generators[i]
            for k in range(j, num_gen):
                # d_ijk = 2 * Tr({T_i, T_j} T_k)
                d_val = 2 * np.trace(anti_comm @ generators[k])
                
                # The result should be real. We take the real part.
                d_val_real = d_val.real
                
                if abs(d_val_real) > TOL:
                    # Round to handle precision issues before adding to set
                    rounded_val = round(d_val_real, 8)
                    d_values.add(rounded_val)

    # Step 3: Output the results
    sorted_values = sorted(list(d_values))
    print(f"Found {len(sorted_values)} distinct non-zero d_ijk values.")
    print("The distinct values are:")
    print(sorted_values)


if __name__ == '__main__':
    # You can change this value to analyze different SU(N) groups.
    # For N=3, the answer is 5.
    # For N=4, the answer is 7.
    # For N=5, the answer is 9.
    # The pattern for 3 <= N <= 6 is 2N-1.
    N_val = 4
    count_distinct_d_values_for_sun(N_val)