import numpy as np

def get_su_n_generators(N):
    """
    Generates a basis for the SU(N) Lie algebra, analogous to the Gell-Mann
    matrices for SU(3). The generators lambda_a are hermitian, traceless, and
    normalized to Tr(lambda_a lambda_b) = 2 * delta_ab.
    These are not the T_a matrices, but are related by T_a = lambda_a / 2.
    """
    if N < 2:
        return []

    num_generators = N**2 - 1
    generators = []

    # Type 1: Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=np.complex128)
            gen[j, k] = 1
            gen[k, j] = 1
            generators.append(gen)

    # Type 2: Antisymmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=np.complex128)
            gen[j, k] = -1j
            gen[k, j] = 1j
            generators.append(gen)

    # Type 3: Diagonal generators
    for l in range(1, N):
        gen = np.zeros((N, N), dtype=np.complex128)
        diag_els = [1.0] * l + [-float(l)] + [0.0] * (N - l - 1)
        # Normalization factor to ensure Tr(gen*gen) = 2
        norm = np.sqrt(2.0 / (l * (l + 1)))
        diag_els_normalized = np.array(diag_els) * norm
        np.fill_diagonal(gen, diag_els_normalized)
        generators.append(gen)

    # Sanity check for the number of generators
    if len(generators) != num_generators:
        raise ValueError("Incorrect number of generators produced.")

    return generators

def count_distinct_d_values(N, precision=8):
    """
    Calculates the number of distinct non-zero values for the d_ijk
    structure constants of SU(N).
    """
    if N < 3:
        # For N=1, there are 0 generators.
        # For N=2, all d_ijk are 0.
        return 0, set()

    # The generators lambda_a are normalized to Tr(lambda_a*lambda_b) = 2*delta_ab
    # T_a = lambda_a / 2, so Tr(T_a*T_b) = 1/2*delta_ab.
    # d_ijk = 4 * Re(Tr(T_i*T_j*T_k)) = 4 * Re(Tr( (lambda_i/2)*(lambda_j/2)*(lambda_k/2) ))
    # d_ijk = 1/2 * Re(Tr(lambda_i*lambda_j*lambda_k))
    lambda_matrices = get_su_n_generators(N)
    num_generators = N**2 - 1

    distinct_values = set()

    # Iterate over all combinations with i <= j <= k
    for i in range(num_generators):
        for j in range(i, num_generators):
            for k in range(j, num_generators):
                # Matrix multiplication
                prod = lambda_matrices[i] @ lambda_matrices[j] @ lambda_matrices[k]
                
                # Calculate the trace
                trace_val = np.trace(prod)

                # d_ijk must be real, so we take the real part.
                d_val = 0.5 * np.real(trace_val)

                # Filter out values that are very close to zero
                if abs(d_val) > 10**(-precision):
                    # Round to handle floating point inaccuracies
                    rounded_val = round(d_val, precision)
                    # Check again after rounding
                    if abs(rounded_val) != 0:
                        distinct_values.add(rounded_val)

    return len(distinct_values), sorted(list(distinct_values))

if __name__ == '__main__':
    # Set the value of N for the SU(N) group
    N = 4

    print(f"Analyzing the structure constants d_ijk for SU({N})...")

    num_distinct_values, values = count_distinct_d_values(N)
    
    print("\nThe distinct non-zero numerical values for d_ijk are:")
    # The user prompt requested to output each number in the final equation.
    # Here we print each distinct value found.
    for val in values:
        print(val)

    print(f"\nFor SU({N}), the number of different numerical values the non-zero d_ijk take is: {num_distinct_values}")
    
    # Final answer in the required format
    # print(f"\n<<<{num_distinct_values}>>>")