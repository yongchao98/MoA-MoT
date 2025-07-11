import numpy as np
import itertools
import sys

def get_su_n_generators(N):
    """
    Generates a set of N^2-1 traceless, Hermitian matrices for SU(N)
    normalized to Tr(L_a L_b) = 2 * delta_ab. These are generalizations
    of the Gell-Mann matrices.

    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        list: A list of N^2-1 numpy arrays representing the generators.
    """
    if N < 2:
        return []

    num_gens = N**2 - 1
    generators = []

    # Off-diagonal symmetric generators
    for k in range(N):
        for l in range(k + 1, N):
            gen = np.zeros((N, N), dtype=np.complex128)
            gen[k, l] = 1
            gen[l, k] = 1
            generators.append(gen)

    # Off-diagonal anti-symmetric generators
    for k in range(N):
        for l in range(k + 1, N):
            gen = np.zeros((N, N), dtype=np.complex128)
            gen[k, l] = -1j
            gen[l, k] = 1j
            generators.append(gen)

    # Diagonal generators
    for m in range(1, N):
        diag_elems = np.zeros(N)
        diag_elems[:m] = 1
        diag_elems[m] = -m
        norm = np.sqrt(2 / (m * (m + 1)))
        gen = norm * np.diag(diag_elems)
        generators.append(gen)

    return generators


def calculate_d_values(N):
    """
    Calculates the unique non-zero values of the symmetric structure 
    constants d_ijk for SU(N).

    Args:
        N (int): The dimension of the special unitary group.

    Returns:
        tuple: A tuple containing the count of unique values and a sorted list 
               of the unique values themselves.
    """
    if not isinstance(N, int) or N < 2:
        raise ValueError("N must be an integer greater than or equal to 2.")

    # For SU(2), all d_ijk are zero.
    if N == 2:
        return 0, []

    # Get the lambda matrices
    generators = get_su_n_generators(N)
    num_gens = len(generators)

    unique_d_values = set()
    # To handle floating point issues, we round to a high precision.
    precision = 10
    tolerance = 10**(-precision + 1)
    
    # Iterate over all combinations of indices i <= j <= k
    indices = range(num_gens)
    for i, j, k in itertools.combinations_with_replacement(indices, 3):
        L_i, L_j, L_k = generators[i], generators[j], generators[k]

        # Calculate d_ijk = (1/4) * Tr({L_i, L_j} L_k)
        anti_comm = L_i @ L_j + L_j @ L_i
        d_val = 0.25 * np.trace(anti_comm @ L_k)

        # The result should be real since the generators are Hermitian.
        d_real = d_val.real
        
        # Add to set if non-zero, after rounding to handle float precision.
        if abs(d_real) > tolerance:
            rounded_val = round(d_real, precision)
            unique_d_values.add(rounded_val)
            
    sorted_values = sorted(list(unique_d_values))
    return len(sorted_values), sorted_values

def main():
    """
    Main function to run the calculation. It prompts for N,
    calculates the d_ijk values, and prints the results.
    """
    try:
        if len(sys.argv) > 1:
            n_val = int(sys.argv[1])
        else:
            n_val_str = input("Enter the value for N (e.g., 3): ")
            n_val = int(n_val_str)

        print(f"\nCalculating for SU({n_val})...")
        count, values = calculate_d_values(n_val)

        print(f"\nFor SU({n_val}), there are {count} different non-zero numerical values for d_ijk.")
        if count > 0:
            print("The unique values are:")
            for val in values:
                # This part is for better symbolic printing of common fractions
                is_printed = False
                for integer in range(1, 5):
                    if np.isclose(val**2, 1/integer):
                        sign = "-" if val < 0 else ""
                        print(f"{val:.6f} \t (Exactly {sign}1/sqrt({integer}))")
                        is_printed = True
                        break
                    for integer2 in range(1, 5):
                        if np.isclose(val**2, integer2/(4*integer)):
                           sign = "-" if val < 0 else ""
                           print(f"{val:.6f} \t (Exactly {sign}sqrt({integer2})/(2*sqrt({integer})))")
                           is_printed = True
                           break
                if not is_printed:
                     print(f"{val:.6f}")


    except ValueError as e:
        print(f"Error: {e}. Please enter a valid integer >= 2.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    # This script requires the numpy library.
    # To run, execute from the command line:
    # python your_script_name.py N
    # where N is an integer (e.g., 3).
    # If N is not provided, it will prompt for input.
    main()
