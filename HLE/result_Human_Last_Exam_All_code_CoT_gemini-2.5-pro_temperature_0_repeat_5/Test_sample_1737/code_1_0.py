import numpy as np

def get_su_n_generators(N):
    """
    Generates a basis of SU(N) generators.
    The basis consists of N^2-1 traceless Hermitian matrices, normalized
    such that Tr(T_a T_b) = 0.5 * delta_ab.
    """
    if N < 2:
        return []

    generators = []
    num_gens = N * N - 1

    # Helper to create a matrix with a 1 at (j,k) and 0 elsewhere
    def E(j, k, n):
        mat = np.zeros((n, n), dtype=complex)
        mat[j, k] = 1
        return mat

    # 1. Symmetric off-diagonal generators
    # N(N-1)/2 of these
    for j in range(N):
        for k in range(j + 1, N):
            # Normalization: Tr( (c(E_jk + E_kj))^2 ) = 0.5
            # c^2 * Tr(E_jj + E_kk) = 0.5 => 2c^2 = 0.5 => c = 0.5
            gen = 0.5 * (E(j, k, N) + E(k, j, N))
            generators.append(gen)

    # 2. Anti-symmetric off-diagonal generators
    # N(N-1)/2 of these
    for j in range(N):
        for k in range(j + 1, N):
            # Normalization: Tr( (c(-i(E_jk - E_kj)))^2 ) = 0.5
            # -c^2 * Tr( -(E_jj + E_kk) ) = 0.5 => 2c^2 = 0.5 => c = 0.5
            gen = -0.5j * (E(j, k, N) - E(k, j, N))
            generators.append(gen)

    # 3. Diagonal generators
    # N-1 of these
    for l in range(1, N):
        # Normalization: Tr( (c * diag(1,..1, -l, 0,..))^2 ) = 0.5
        # c^2 * (l*1^2 + (-l)^2) = 0.5 => c^2 * l(l+1) = 0.5
        # c = 1 / sqrt(2*l*(l+1))
        norm = 1.0 / np.sqrt(2 * l * (l + 1))
        gen = np.zeros((N, N), dtype=complex)
        for m in range(l):
            gen[m, m] = norm
        gen[l, l] = -l * norm
        generators.append(gen)

    return generators

def solve_for_n(N):
    """
    Calculates the number of unique non-zero d_ijk values for SU(N).
    """
    print(f"Calculating for SU({N})...")

    if N < 2:
        print("N must be >= 2.")
        return
        
    if N == 2:
        print("For SU(2), all d_ijk coefficients are zero.")
        print("Number of different non-zero values: 0")
        print("<<<0>>>")
        return

    generators = get_su_n_generators(N)
    num_gens = len(generators)
    
    d_values = set()
    
    # Iterate over all combinations i <= j <= k
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                Ti = generators[i]
                Tj = generators[j]
                Tk = generators[k]
                
                # Anti-commutator {Ti, Tj}
                anti_comm = Ti @ Tj + Tj @ Ti
                
                # d_ijk = 2 * Tr({Ti, Tj} Tk)
                d_ijk = 2 * np.trace(anti_comm @ Tk)
                
                # Check if the value is non-zero within a tolerance
                if abs(d_ijk.real) > 1e-9:
                    # Round to handle floating point inaccuracies
                    rounded_val = round(d_ijk.real, 8)
                    d_values.add(rounded_val)

    print(f"Found {len(d_values)} different non-zero d_ijk values for SU({N}).")
    
    # The final equation part: printing each number
    if d_values:
        sorted_values = sorted(list(d_values))
        print("The unique non-zero values are:")
        # Create the equation string
        equation_str = "{" + ", ".join(map(str, sorted_values)) + "}"
        print(equation_str)
    
    # The analytical result for N>=3 is 4*N - 7
    analytical_result = 4 * N - 7
    print(f"\nAnalytical formula for N>=3 gives: 4*N - 7 = 4*{N} - 7 = {analytical_result}")
    
    # Final answer format
    print(f"<<<{len(d_values)}>>>")


if __name__ == '__main__':
    # Set the value of N to solve for
    # For example, N=3, N=4, etc.
    N = 3
    solve_for_n(N)
