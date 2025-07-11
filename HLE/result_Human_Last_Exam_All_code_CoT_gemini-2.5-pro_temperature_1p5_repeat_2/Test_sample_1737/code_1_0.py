import numpy as np

def get_gell_mann_matrices(n):
    """
    Generates the N^2 - 1 generalized Gell-Mann matrices for SU(N).
    These matrices (lambda_a) are a basis of the N x N traceless Hermitian matrices.
    They are normalized such that Tr(lambda_a * lambda_b) = 2 * delta_ab.
    """
    if not isinstance(n, int) or n < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    num_matrices = n**2 - 1
    matrices = []

    # Type 1: Symmetric off-diagonal matrices
    # There are N*(N-1)/2 of these.
    for i in range(n):
        for j in range(i + 1, n):
            matrix = np.zeros((n, n), dtype=complex)
            matrix[i, j] = 1
            matrix[j, i] = 1
            matrices.append(matrix)

    # Type 2: Antisymmetric off-diagonal matrices
    # There are N*(N-1)/2 of these.
    for i in range(n):
        for j in range(i + 1, n):
            matrix = np.zeros((n, n), dtype=complex)
            matrix[i, j] = -1j
            matrix[j, i] = 1j
            matrices.append(matrix)

    # Type 3: Diagonal matrices
    # There are N-1 of these.
    for k in range(1, n):
        matrix = np.zeros((n, n), dtype=complex)
        val = np.sqrt(2.0 / (k * (k + 1)))
        for i in range(k):
            matrix[i, i] = val
        matrix[k, k] = -k * val
        matrices.append(matrix)

    return matrices

def calculate_d_values(n):
    """
    Calculates the set of unique non-zero d_ijk values for SU(N).
    """
    try:
        lambdas = get_gell_mann_matrices(n)
    except ValueError as e:
        print(f"Error: {e}")
        return

    num_gens = n**2 - 1
    d_values = set()
    
    # We only need to iterate over i <= j <= k due to symmetry
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                l_i, l_j, l_k = lambdas[i], lambdas[j], lambdas[k]
                
                # Calculate anticommutator {l_i, l_j}
                anticomm = np.dot(l_i, l_j) + np.dot(l_j, l_i)
                
                # Calculate d_ijk = (1/4) * Tr({l_i, l_j} * l_k)
                d_val = 0.25 * np.trace(np.dot(anticomm, l_k))
                
                # d_ijk must be real. Take the real part to discard numerical noise.
                d_val_real = np.real(d_val)
                
                # Add to set if non-zero (within a tolerance)
                if abs(d_val_real) > 1e-9:
                    d_values.add(round(d_val_real, 8))
                    
    return sorted(list(d_values))

if __name__ == '__main__':
    # Set the value of N for SU(N)
    N = 4 

    if N < 2:
         print(f"For SU({N}), the concept of d_ijk is not standard.")
    elif N == 2:
        # For SU(2), all d_ijk are zero.
        num_distinct_values = 0
        distinct_values = []
        print(f"For SU({N}), all d_ijk constants are zero.")
        print(f"Number of distinct non-zero values: {num_distinct_values}")

    else:
        distinct_values = calculate_d_values(N)
        num_distinct_values = len(distinct_values)

        print(f"For SU({N}), the number of distinct non-zero d_ijk values is: {num_distinct_values}")
        print("The distinct numerical values are:")
        for val in distinct_values:
            print(f"{val:.8f}")
