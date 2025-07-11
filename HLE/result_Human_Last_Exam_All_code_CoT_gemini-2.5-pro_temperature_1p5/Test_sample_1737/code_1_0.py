import numpy as np

def get_su_n_generators(N):
    """
    Generates a basis of SU(N) generators T_a, normalized as Tr(T_a T_b) = 0.5 * delta_ab.
    These are (N^2 - 1) traceless Hermitian matrices of size N x N.
    """
    if N < 2:
        return []

    generators = []
    
    # Off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            # Symmetric generators (real entries)
            s_mat = np.zeros((N, N), dtype=complex)
            s_mat[j, k] = 0.5
            s_mat[k, j] = 0.5
            generators.append(s_mat)
            
            # Anti-symmetric generators (imaginary entries)
            a_mat = np.zeros((N, N), dtype=complex)
            a_mat[j, k] = -0.5j
            a_mat[k, j] = 0.5j
            generators.append(a_mat)

    # Diagonal generators
    for l in range(1, N):
        d_mat = np.zeros((N, N), dtype=complex)
        # Normalization coefficient
        coeff = 1.0 / np.sqrt(2 * l * (l + 1))
        for i in range(l):
            d_mat[i, i] = coeff
        d_mat[l, l] = -l * coeff
        generators.append(d_mat)
        
    return generators

def solve_su_n_d_constants():
    """
    Calculates the number of unique non-zero values for the symmetric structure constants
    d_ijk of SU(N) for a given N, and prints the result.
    """
    try:
        N_str = input("Enter the value of N for SU(N): ")
        N = int(N_str)
        if N < 2:
            print("N must be an integer greater than or equal to 2.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    print(f"\nStarting the calculation for SU({N})...")
    
    # 1. Generate the SU(N) generators
    generators = get_su_n_generators(N)
    dim = len(generators)
    
    # Handle the trivial case for SU(2), where all d_ijk are zero.
    if N == 2:
        num_unique_values = 0
        unique_d_values = set()
    else:
        # 2. Iterate and calculate d_ijk values
        unique_d_values = set()
        TOLERANCE = 1e-9 # To handle floating point comparisons
        
        # We use i <= j <= k due to the total symmetry of d_ijk
        for i in range(dim):
            for j in range(i, dim):
                for k in range(j, dim):
                    Ti, Tj, Tk = generators[i], generators[j], generators[k]
                    
                    # Calculate d_ijk = 2 * Tr(T_i {T_j, T_k})
                    # where {T_j, T_k} = Tj @ Tk + Tk @ Tj
                    anticommutator = np.dot(Tj, Tk) + np.dot(Tk, Tj)
                    d_ijk_complex = 2 * np.trace(np.dot(Ti, anticommutator))
                    
                    # The result should be real; take the real part to discard numerical noise.
                    d_ijk = np.real(d_ijk_complex)
                    
                    # 3. Store unique non-zero values
                    if abs(d_ijk) > TOLERANCE:
                        # Round to a fixed precision to group similar floating point numbers
                        rounded_val = round(d_ijk, 8)
                        unique_d_values.add(rounded_val)

        num_unique_values = len(unique_d_values)
        
    # 4. Print the result
    print(f"\nFor SU({N}), the totally symmetric structure constants d_ijk take {num_unique_values} different non-zero numerical values.")
    
    if num_unique_values > 0:
        sorted_values = sorted(list(unique_d_values))
        print("These values are:")
        print(sorted_values)

# Execute the main function to solve the problem
solve_su_n_d_constants()
