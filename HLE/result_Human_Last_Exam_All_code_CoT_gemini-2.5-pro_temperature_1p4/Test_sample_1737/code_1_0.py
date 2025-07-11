import numpy as np
import itertools

def get_gell_mann_matrices(N):
    """
    Generates the N^2 - 1 generalized Gell-Mann matrices for SU(N).
    The matrices follow the standard ordering convention.
    """
    if N < 2:
        return []
    
    d = N**2 - 1
    mats = [np.zeros((N, N), dtype=complex) for _ in range(d)]

    # 1. Identify indices for diagonal matrices
    # For SU(N), the diagonal lambda matrices are at indices m^2-2 for m=2,...,N (0-based)
    diag_indices = {m**2 - 2 for m in range(2, N + 1)}

    # 2. Generate and place off-diagonal matrices
    current_mat_idx = 0
    for j in range(N):
        for k in range(j + 1, N):
            # Place symmetric matrix
            while current_mat_idx in diag_indices:
                current_mat_idx += 1
            s = np.zeros((N, N), dtype=complex)
            s[j, k] = 1
            s[k, j] = 1
            mats[current_mat_idx] = s
            current_mat_idx += 1

            # Place antisymmetric matrix
            while current_mat_idx in diag_indices:
                current_mat_idx += 1
            a = np.zeros((N, N), dtype=complex)
            a[j, k] = -1j
            a[k, j] = 1j
            mats[current_mat_idx] = a
            current_mat_idx += 1

    # 3. Generate and place diagonal matrices
    for m in range(2, N + 1):
        idx = m**2 - 2
        mat = np.zeros((N, N), dtype=complex)
        # Normalization factor sqrt(2/(m(m-1)))
        factor = np.sqrt(2.0 / (m * (m - 1)))
        for i in range(m - 1):
            mat[i, i] = factor
        mat[m - 1, m - 1] = -(m - 1) * factor
        mats[idx] = mat
            
    return mats

def solve_for_su_n(N):
    """
    Calculates the number of distinct non-zero values for the symmetric
    structure constants d_ijk of SU(N).
    """
    if N < 2:
        print(f"SU({N}) is trivial or not well-defined in this context.")
        print("Number of different non-zero values: 0")
        return 0

    print(f"Analyzing for SU({N})...")
    
    # Step 1: Get the generator matrices (lambda matrices)
    lambdas = get_gell_mann_matrices(N)
    d = N**2 - 1

    # Step 2: Iterate through combinations and calculate d_ijk
    # We use a set to store unique values
    unique_d_values = set()
    
    # To be efficient, we only loop over i <= j <= k due to symmetry
    indices = range(d)
    for i, j, k in itertools.combinations_with_replacement(indices, 3):
        # Formula: d_ijk = (1/4) * Tr({lambda_i, lambda_j} * lambda_k)
        lambda_i, lambda_j, lambda_k = lambdas[i], lambdas[j], lambdas[k]
        
        # {lambda_i, lambda_j} = lambda_i * lambda_j + lambda_j * lambda_i
        anti_commutator = np.dot(lambda_i, lambda_j) + np.dot(lambda_j, lambda_i)
        
        # Calculate d_ijk
        d_ijk_val = 0.25 * np.trace(np.dot(anti_commutator, lambda_k))
        
        # The result should be real. Take the real part to discard small imaginary noise.
        d_ijk_val = np.real(d_ijk_val)

        # Step 3: Add to the set if non-zero (using a tolerance for floating point)
        if not np.isclose(d_ijk_val, 0, atol=1e-9):
            # Round the value to a standard precision to group similar float results
            rounded_val = round(d_ijk_val, 10)
            unique_d_values.add(rounded_val)

    # Step 4: Print the results
    print("The different numerical values of the non-zero d_ijk are:")
    sorted_values = sorted(list(unique_d_values))
    for val in sorted_values:
        print(val)
    
    num_unique_values = len(unique_d_values)
    print(f"\nFor SU({N}), the number of different non-zero d_ijk values is: {num_unique_values}")
    return num_unique_values

if __name__ == '__main__':
    # Set the value of N for the SU(N) group you want to analyze.
    # For example, N=2, N=3, N=4, etc.
    # The calculation can become slow for N > 5.
    N = 3
    
    solve_for_su_n(N)
    
    # As an example, let's also run for N=4
    # print("\n" + "="*40 + "\n")
    # solve_for_su_n(4)

