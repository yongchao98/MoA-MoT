import numpy as np

def get_su_n_generators(N):
    """
    Generates a basis for the SU(N) Lie algebra, known as the
    generalized Gell-Mann matrices. They are traceless, Hermitian,
    and satisfy Tr(L_a L_b) = 2 * delta_ab.
    
    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        list: A list of NxN numpy arrays representing the generators.
    """
    if not isinstance(N, int) or N < 2:
        return []

    num_gens = N**2 - 1
    generators = []

    # Type 1: Symmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = 1
            mat[k, j] = 1
            generators.append(mat)

    # Type 2: Antisymmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            mat = np.zeros((N, N), dtype=complex)
            mat[j, k] = -1j
            mat[k, j] = 1j
            generators.append(mat)

    # Type 3: Diagonal matrices
    for l in range(1, N):
        mat = np.zeros((N, N), dtype=complex)
        for m in range(l):
            mat[m, m] = 1
        mat[l, l] = -l
        mat *= np.sqrt(2.0 / (l * (l + 1)))
        generators.append(mat)

    # Verification that we have the correct number of generators
    assert len(generators) == num_gens, "Incorrect number of generators created."
    
    return generators

def calculate_d_values(N):
    """
    Calculates the distinct non-zero d_ijk values for SU(N).
    
    Args:
        N (int): The dimension of the special unitary group.
    """
    if N < 3:
        print(f"\nFor SU({N}), all d_ijk coefficients are zero.")
        print("Number of distinct non-zero values: 0")
        return

    generators = get_su_n_generators(N)
    num_gens = len(generators)
    
    d_values = set()
    # To handle floating point inaccuracies, we will round the results
    precision = 8

    # d_ijk is totally symmetric, so we iterate with i <= j <= k
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                # Calculate anticommutator {L_i, L_j}
                anticomm = np.dot(generators[i], generators[j]) + np.dot(generators[j], generators[i])
                
                # Calculate d_ijk = (1/4) * Tr({L_i, L_j} L_k)
                d_ijk = 0.25 * np.trace(np.dot(anticomm, generators[k]))
                
                # The result must be real; np.real extracts it and discards tiny imaginary errors
                d_ijk_real = np.real(d_ijk)

                # Add non-zero values to a set to count unique ones
                if abs(d_ijk_real) > 10**(-precision):
                    d_values.add(round(d_ijk_real, precision))

    num_distinct_values = len(d_values)
    print(f"\nFor SU({N}), the totally symmetric structure constants d_ijk take {num_distinct_values} distinct non-zero numerical values.")
    
    if num_distinct_values > 0:
        print("These values are:")
        sorted_values = sorted(list(d_values))
        for val in sorted_values:
            print(val)

def main():
    """
    Main function to get user input and run the calculation.
    """
    try:
        n_input_str = input("Enter the value for N (must be >= 2): ")
        N = int(n_input_str)
        if N < 2:
            print("Invalid input. N must be an integer greater than or equal to 2.")
        else:
            calculate_d_values(N)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")

if __name__ == "__main__":
    main()
