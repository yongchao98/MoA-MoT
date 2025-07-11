import numpy as np
from itertools import combinations_with_replacement

def get_su_n_generators(N):
    """
    Generates the N^2 - 1 traceless Hermitian generators of SU(N)
    normalized such that Tr(T_a T_b) = (1/2) * delta_ab.
    
    Args:
        N (int): The dimension of the group SU(N).

    Returns:
        list: A list of numpy arrays representing the generators.
    """
    if not isinstance(N, int) or N < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    generators = []
    
    # Type 1: Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            # Normalization factor of 1/2 from Tr(T^2) = 1/2
            gen[j, k] = 0.5
            gen[k, j] = 0.5
            generators.append(gen)
    
    # Type 2: Anti-symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            gen = np.zeros((N, N), dtype=complex)
            # Normalization factor of 1/2 from Tr(T^2) = 1/2
            gen[j, k] = -0.5j
            gen[k, j] = 0.5j
            generators.append(gen)

    # Type 3: Diagonal generators
    for l in range(1, N):
        diag_vals = np.zeros(N, dtype=float)
        diag_vals[:l] = 1.0
        diag_vals[l] = -float(l)
        
        # Normalization constant C from Tr((C*diag_vals)^2) = 1/2
        # C^2 * (l*1^2 + (-l)^2) = 1/2 => C = 1 / sqrt(2*l*(l+1))
        norm = np.sqrt(2 * l * (l + 1))
        gen = np.diag(diag_vals / norm)
        generators.append(gen)
        
    return generators

def calculate_and_count_d_ijk(N):
    """
    Calculates the d_ijk structure constants for SU(N) and counts the
    number of unique non-zero values.
    
    Args:
        N (int): The dimension of the group SU(N).
    """
    try:
        T = get_su_n_generators(N)
        num_gens = len(T)
        
        d_values = set()
        
        # d_ijk is totally symmetric, so iterate over i <= j <= k
        indices = range(num_gens)
        for i, j, k in combinations_with_replacement(indices, 3):
            T_i, T_j, T_k = T[i], T[j], T[k]
            
            # Formula: d_ijk = 2 * Tr(T_i T_j T_k + T_j T_i T_k)
            term1 = T_i @ T_j @ T_k
            term2 = T_j @ T_i @ T_k
            
            # d_ijk must be real. Use np.real to discard tiny imaginary parts from float errors.
            d_val = np.real(2 * np.trace(term1 + term2))
            
            # Consider a value non-zero if its magnitude is greater than a small threshold
            if abs(d_val) > 1e-9:
                # Round to a fixed precision to group numerically close values
                d_values.add(round(d_val, 8))
        
        print(f"For SU({N}), the unique non-zero d_ijk values are:")
        if not d_values:
            print("None")
        else:
            # Sort the values for consistent output
            sorted_values = sorted(list(d_values))
            print("{", end="")
            print(*sorted_values, sep=", ", end="")
            print("}")
        
        print(f"\nThe number of different numerical values the non-zero d_ijk take is:")
        print(len(d_values))
        return len(d_values)

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    try:
        N_str = input("Enter the value of N for SU(N) (e.g., 2, 3, 4): ")
        N = int(N_str)
        result = calculate_and_count_d_ijk(N)
        if result is not None:
            print(f"<<<{result}>>>")
    except ValueError:
        print("Invalid input. Please enter an integer.")
