import numpy as np

def count_d_values(N):
    """
    Calculates the number of unique non-zero d_ijk values for SU(N).

    Args:
        N (int): The dimension of the special unitary group SU(N).

    Returns:
        int: The number of unique non-zero d_ijk values.
    """
    if N < 2:
        print("N must be an integer greater than or equal to 2.")
        return 0
        
    if N == 2:
        # For SU(2), all d_ijk are zero.
        print("For SU(2), all d_ijk coefficients are zero.")
        print("Number of different non-zero values: 0")
        return 0

    dim = N * N - 1
    generators = []

    # 1. Generate the SU(N) generators (T_a = lambda_a / 2)

    # Off-diagonal generators
    # Symmetric (real)
    for j in range(N):
        for k in range(j + 1, N):
            # lambda_jk^S has 1 at (j,k) and (k,j)
            # Tr((lambda_jk^S)^2) = 2
            l_mat = np.zeros((N, N), dtype=complex)
            l_mat[j, k] = 1.0
            l_mat[k, j] = 1.0
            generators.append(l_mat / 2.0)

    # Anti-symmetric (imaginary)
    for j in range(N):
        for k in range(j + 1, N):
            # lambda_jk^A has -i at (j,k) and i at (k,j)
            # Tr((lambda_jk^A)^2) = 2
            l_mat = np.zeros((N, N), dtype=complex)
            l_mat[j, k] = -1.0j
            l_mat[k, j] = 1.0j
            generators.append(l_mat / 2.0)

    # Diagonal generators
    for l in range(1, N):
        # lambda_l^D = sqrt(2/(l(l+1))) * diag(1, ..., 1, -l, 0, ...)
        # Tr((lambda_l^D)^2) = 2
        l_mat = np.zeros((N, N), dtype=complex)
        norm = np.sqrt(2.0 / (l * (l + 1.0)))
        for i in range(l):
            l_mat[i, i] = norm
        l_mat[l, l] = -l * norm
        generators.append(l_mat / 2.0)

    # 2. Calculate d_ijk values and find unique ones
    unique_d_values = set()
    epsilon = 1e-9

    for i in range(dim):
        for j in range(i, dim):
            for k in range(j, dim):
                T_i = generators[i]
                T_j = generators[j]
                T_k = generators[k]

                # d_ijk = 2 * Tr({T_i, T_j} T_k) = 2 * Tr((T_i T_j + T_j T_i) T_k)
                d_val = 2 * np.trace((T_i @ T_j + T_j @ T_i) @ T_k)

                # The result should be real. We take the real part.
                d_val_real = d_val.real

                if abs(d_val_real) > epsilon:
                    # Round to handle floating point inaccuracies
                    rounded_val = round(d_val_real, 7)
                    unique_d_values.add(rounded_val)

    # 3. Print the results
    count = len(unique_d_values)
    print(f"\nFor SU(N={N}), the calculation is complete.")
    print(f"The number of different non-zero d_ijk values is: {count}")
    
    if count > 0:
        print("\nThe unique values found are:")
        sorted_values = sorted(list(unique_d_values))
        for val in sorted_values:
            print(f"{val: .8f}")
            
    return count

if __name__ == '__main__':
    try:
        n_input = input("Enter the value of N for SU(N): ")
        N_val = int(n_input)
        count_d_values(N_val)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
