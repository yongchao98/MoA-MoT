import numpy as np
from itertools import combinations_with_replacement

def get_sun_generators(N):
    """
    Generates a basis of SU(N) generators (generalized Gell-Mann matrices).
    The generators T_a are normalized such that Tr(T_a T_b) = (1/2) * delta_ab.
    
    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        list: A list of N^2 - 1 generator matrices.
    """
    if not isinstance(N, int) or N < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    generators = []

    # Off-diagonal generators
    for p in range(N):
        for q in range(p + 1, N):
            # Symmetric generators
            S_pq = np.zeros((N, N), dtype=complex)
            S_pq[p, q] = 1
            S_pq[q, p] = 1
            generators.append(S_pq / 2.0)
            
            # Anti-symmetric generators
            A_pq = np.zeros((N, N), dtype=complex)
            A_pq[p, q] = -1j
            A_pq[q, p] = 1j
            generators.append(A_pq / 2.0)

    # Diagonal generators
    for l in range(1, N):
        diag_vals = np.zeros(N, dtype=float)
        diag_vals[:l] = 1
        diag_vals[l] = -l
        norm = np.sqrt(2 * l * (l + 1))
        D_l = np.diag(diag_vals) / norm
        generators.append(D_l)
        
    return generators

def find_d_values(N):
    """
    Calculates the set of unique non-zero d_ijk values for SU(N).
    
    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        list: A sorted list of unique non-zero d_ijk values.
    """
    if N < 2:
        return []
    # For SU(2), all d_ijk are zero.
    if N == 2:
        return []

    generators = get_sun_generators(N)
    dim = len(generators)
    d_values = set()
    
    # Iterate through all combinations of indices i, j, k with i <= j <= k
    indices = range(dim)
    for i, j, k in combinations_with_replacement(indices, 3):
        Ti = generators[i]
        Tj = generators[j]
        Tk = generators[k]
        
        # d_ijk = 2 * Tr({Ti, Tj} Tk) = 2 * Tr(Ti Tj Tk + Tj Ti Tk)
        d_ijk_val = 2 * np.trace(Ti @ Tj @ Tk + Tj @ Ti @ Tk)
        
        # The result must be real since d_ijk are real.
        d_ijk_val = np.real(d_ijk_val)
        
        # Check if the value is non-zero within a tolerance
        if abs(d_ijk_val) > 1e-9:
            # Round to handle floating point inaccuracies before adding to set
            rounded_val = np.round(d_ijk_val, 8)
            d_values.add(rounded_val)
            
    return sorted(list(d_values))

def main():
    """
    Main function to run the calculation and print results.
    """
    try:
        N_str = input("Enter the value of N for SU(N) (e.g., 3, 4, 5): ")
        N = int(N_str)
        
        print(f"\nCalculating for SU({N})...")
        unique_values = find_d_values(N)
        num_values = len(unique_values)
        
        print(f"\nFor SU({N}), there are {num_values} different numerical values for the non-zero d_ijk constants.")
        
        if num_values > 0:
            print("The unique values are:")
            # To fulfill the instruction "output each number in the final equation",
            # we print each of the unique values found.
            for val in unique_values:
                print(f"{val: .8f}")

    except ValueError as e:
        print(f"Error: {e}. Please enter an integer N >= 2.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()