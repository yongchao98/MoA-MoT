import numpy as np

def get_sun_generators(N):
    """
    Constructs a set of generators for SU(N) normalized to Tr(T_a T_b) = 0.5 * delta_ab.
    These are the generalized Gell-Mann matrices, divided by 2.
    """
    if N < 2:
        return []
    
    generators = []
    
    # Type 1: Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            # T = c * (E_jk + E_kj). Tr(T^2) = 2c^2. For Tr(T^2)=1/2, c=1/2.
            mat[j, k] = 0.5
            mat[k, j] = 0.5
            generators.append(mat)
    
    # Type 2: Antisymmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            # T is Hermitian, T = c*i*(E_kj - E_jk). Tr(T^2)=2c^2. For Tr(T^2)=1/2, c=1/2.
            mat[j, k] = -0.5j
            mat[k, j] = 0.5j
            generators.append(mat)

    # Type 3: Diagonal generators
    for l in range(1, N):
        mat = np.zeros((N, N), dtype=complex)
        # T is diagonal and real.
        # diag(1, ..., 1 (l terms), -l, 0, ...), norm c=1/sqrt(2*l*(l+1)) for Tr(T^2)=1/2
        diag_vec = np.zeros(N, dtype=float)
        diag_vec[0:l] = 1
        diag_vec[l] = -l
        norm = np.sqrt(2 * l * (l + 1))
        diag_vec /= norm
        np.fill_diagonal(mat, diag_vec)
        generators.append(mat)
    
    return generators

def calculate_d_values(N):
    """
    Calculates the number of unique, non-zero d_ijk values for SU(N).
    """
    try:
        N = int(N)
        if N < 2:
            print("N must be an integer greater than or equal to 2.")
            return
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer N >= 2.")
        return

    generators = get_sun_generators(N)
    num_gen = len(generators)
    
    d_values = set()
    
    # A small tolerance for checking if a value is non-zero
    TOLERANCE = 1e-9
    # Round to a fixed number of decimal places to group floating point numbers
    PRECISION = 8

    # Iterate over all combinations of i, j, k with i <= j <= k
    for i in range(num_gen):
        for j in range(i, num_gen):
            for k in range(j, num_gen):
                Ti, Tj, Tk = generators[i], generators[j], generators[k]
                
                # d_ijk = 2 * Tr( (Ti Tj + Tj Ti) Tk )
                d_val_complex = 2 * np.trace( (Ti @ Tj + Tj @ Ti) @ Tk )
                
                # The d_ijk must be real. The imaginary part is due to numerical error.
                d_val = np.real(d_val_complex)
                
                # If the value is non-zero within the tolerance, add it to the set.
                if abs(d_val) > TOLERANCE:
                    d_values.add(round(d_val, PRECISION))
    
    # Final Output
    print(f"For SU({N}), the totally symmetric structure constants d_ijk can take the following unique non-zero values:")
    sorted_values = sorted(list(d_values))
    print(sorted_values)
    
    num_distinct_values = len(d_values)
    print(f"\nThe number of different numerical values for the non-zero d_ijk is:")
    print(num_distinct_values)
    
    return num_distinct_values

if __name__ == '__main__':
    try:
        n_input = input("Enter the value for N (e.g., 3): ")
        result = calculate_d_values(n_input)
        if result is not None:
            # The final answer is wrapped according to the format specification
            print(f"<<<{result}>>>")
    except Exception as e:
        print(f"An error occurred: {e}")