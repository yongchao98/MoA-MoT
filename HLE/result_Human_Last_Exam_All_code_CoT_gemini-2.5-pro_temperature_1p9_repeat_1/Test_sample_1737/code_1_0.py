import numpy as np

def solve_su_n_d_values():
    """
    Calculates the number of distinct non-zero values for the totally symmetric
    structure constants d_ijk of SU(N) for a user-provided N.
    """
    try:
        n_str = input("Enter the value for N (integer >= 2): ")
        N = int(n_str)
        if N < 2:
            print("N must be an integer greater than or equal to 2.")
            return
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
        return

    # Step 1: Generate the SU(N) generators, normalized to Tr(T_a T_b) = 0.5 * delta_ab
    generators = []
    
    # Pre-build basis matrices E_ij (1 at (i,j), 0 elsewhere) for convenience
    E = {}
    for j in range(N):
        for k in range(N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = 1
            E[(j, k)] = mat

    # Type 1: Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = 0.5 * (E[(j, k)] + E[(k, j)])
            generators.append(gen)

    # Type 2: Antisymmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = -0.5j * (E[(j, k)] - E[(k, j)])
            generators.append(gen)

    # Type 3: Diagonal generators
    for l in range(1, N):
        diag_entries = np.zeros(N, dtype=float)
        diag_entries[:l] = 1
        diag_entries[l] = -l
        norm_factor = 1.0 / np.sqrt(2 * l * (l + 1))
        gen = np.diag(diag_entries) * norm_factor
        generators.append(gen)

    num_gens = len(generators)
    
    # Step 2: Calculate d_ijk = 4 * Re(Tr(T_i T_j T_k)) for all i <= j <= k
    d_values = set()
    # Precision for rounding to handle floating point inaccuracies
    precision = 10 
    
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                # Using np.linalg.multi_dot is efficient for product of multiple matrices
                trace_val = np.trace(np.linalg.multi_dot([generators[i], generators[j], generators[k]]))
                
                # The d_ijk coefficient is 4 * Re(trace)
                d_val = 4 * np.real(trace_val)

                # Check if the value is non-zero within a tolerance
                if not np.isclose(d_val, 0, atol=1e-9):
                    d_val_rounded = round(d_val, precision)
                    d_values.add(d_val_rounded)

    # Step 3: Print the results
    num_unique_values = len(d_values)
    print(f"\nFor SU({N}), there are {num_unique_values} distinct non-zero numerical values for d_ijk.")
    
    if num_unique_values > 0:
        print("The distinct values are:")
        # Sort for consistent ordering
        sorted_values = sorted(list(d_values))
        for val in sorted_values:
            print(val)

# Execute the main function
solve_su_n_d_values()