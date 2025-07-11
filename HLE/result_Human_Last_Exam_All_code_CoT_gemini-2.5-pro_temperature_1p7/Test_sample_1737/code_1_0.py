import numpy as np
import itertools

def solve_su_n_d_values(N):
    """
    Calculates the number of unique non-zero values for the symmetric
    structure constants d_ijk of SU(N).

    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        None. Prints the results.
    """
    if not isinstance(N, int) or N < 2:
        print("N must be an integer greater than or equal to 2.")
        return

    print(f"Calculating for SU({N})...")

    # 1. Generate the N^2 - 1 generators (generalized Gell-Mann matrices)
    # The generators T_a = lambda_a / 2 are normalized such that Tr(T_a T_b) = 0.5 * delta_ab.
    
    lambdas = []
    
    # Type 1: Symmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            gen[j, k] = 1.0
            gen[k, j] = 1.0
            lambdas.append(gen)

    # Type 2: Antisymmetric off-diagonal matrices
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            gen[j, k] = -1.0j
            gen[k, j] = 1.0j
            lambdas.append(gen)

    # Type 3: Diagonal matrices
    for l in range(1, N):
        gen = np.zeros((N, N), dtype=complex)
        norm_factor = np.sqrt(2.0 / (l * (l + 1)))
        for m in range(l):
            gen[m, m] = norm_factor
        gen[l, l] = -l * norm_factor
        lambdas.append(gen)

    generators = [lam / 2.0 for lam in lambdas]
    num_gen = len(generators)

    # For SU(2), all d_ijk are 0
    if N == 2:
        print("For SU(2), all d_ijk constants are zero.")
        print("Number of different non-zero values: 0")
        print("The set of unique non-zero values is: {}")
        # Final answer format for this specific case
        print("\n<<<0>>>")
        return

    # 2. Calculate d_ijk values and find the unique non-zero ones
    d_values = set()
    
    # Exploit symmetry by iterating with i <= j <= k
    indices = range(num_gen)
    for i, j, k in itertools.combinations_with_replacement(indices, 3):
        
        # d_ijk = Tr({T_i, T_j} T_k)
        anti_commutator = np.dot(generators[i], generators[j]) + np.dot(generators[j], generators[i])
        d_val = np.trace(np.dot(anti_commutator, generators[k]))

        # The result must be real. Take real part to handle floating point inaccuracies.
        d_val = d_val.real
        
        # Check if the value is non-zero within a tolerance
        if not np.isclose(d_val, 0):
            # Round to group numerically identical values
            d_values.add(round(d_val, 10))

    # 3. Print the results
    sorted_values = sorted(list(d_values))
    count = len(sorted_values)

    print(f"\nFor SU({N}), the number of different non-zero numerical values for d_ijk is: {count}")
    print("The unique values are:")
    for val in sorted_values:
        print(val)
    
    # For the final answer submission as requested by the prompt format
    # we return the value for the chosen N=3.
    if N == 3:
        print(f"\n<<<{count}>>>")


if __name__ == '__main__':
    # To solve for a different SU(N), change the value of N below.
    # The default is N=3, the first non-trivial case.
    N_val = 3
    solve_su_n_d_values(N_val)