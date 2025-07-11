import numpy as np
from itertools import combinations_with_replacement

def solve_for_su_n():
    """
    Calculates the number of different numerical values for the non-zero symmetric
    structure constants d_ijk of SU(N).
    """
    try:
        N_input = input("Enter the value for N (e.g., 3 for SU(3)): ")
        N = int(N_input)
        if N < 2:
            print("N must be an integer greater than or equal to 2.")
            return
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer.")
        return

    print(f"\nCalculating for SU({N})...")

    # SU(2) is a special case where all d_ijk are zero.
    if N == 2:
        print("For SU(2), all d_ijk constants are zero.")
        num_unique_values = 0
    else:
        # 1. Construct SU(N) generators, normalized as Tr(T_a T_b) = 0.5 * delta_ab
        generators = []
        
        # A. Off-diagonal symmetric matrices (N*(N-1)/2 of them)
        for j in range(N):
            for k in range(j + 1, N):
                T = np.zeros((N, N), dtype=complex)
                T[j, k] = 1
                T[k, j] = 1
                generators.append(T / 2.0)

        # B. Off-diagonal anti-symmetric matrices (N*(N-1)/2 of them)
        for j in range(N):
            for k in range(j + 1, N):
                T = np.zeros((N, N), dtype=complex)
                T[j, k] = -1j
                T[k, j] = 1j
                generators.append(T / 2.0)

        # C. Diagonal matrices (N-1 of them)
        for l in range(1, N):
            diag_entries = [1.0] * l + [-float(l)] + [0.0] * (N - l - 1)
            T = np.diag(diag_entries).astype(complex)
            norm = np.sqrt(2 * l * (l + 1))
            generators.append(T / norm)
        
        num_gen = len(generators)
        unique_d_values = set()
        TOLERANCE = 1e-9

        indices = range(num_gen)
        
        # 2. Iterate over all unique combinations (i,j,k) due to symmetry
        for i, j, k in combinations_with_replacement(indices, 3):
            Ti = generators[i]
            Tj = generators[j]
            Tk = generators[k]
            
            # 3. Calculate d_ijk = 2 * Tr({Ti, Tj} * Tk)
            anticommutator = np.dot(Ti, Tj) + np.dot(Tj, Ti)
            d_val = np.real(2 * np.trace(np.dot(anticommutator, Tk)))

            # 4. Store unique non-zero values, handling float precision
            if abs(d_val) > TOLERANCE:
                unique_d_values.add(round(d_val, 8))

        num_unique_values = len(unique_d_values)

    print("-" * 30)
    print(f"Result for SU({N}):")
    print(f"The number of different numerical values for non-zero d_ijk is: {num_unique_values}")
    print("-" * 30)

if __name__ == '__main__':
    solve_for_su_n()
