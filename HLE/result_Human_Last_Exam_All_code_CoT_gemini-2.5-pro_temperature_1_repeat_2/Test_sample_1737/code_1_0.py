import numpy as np
import sys

def generate_gell_mann(N):
    """
    Generates the N^2 - 1 generalized Gell-Mann matrices for SU(N).
    These are a standard basis of traceless Hermitian matrices.
    """
    if N < 2:
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
        norm = np.sqrt(2 / (l * (l + 1)))
        for m in range(l):
            mat[m, m] = 1
        mat[l, l] = -l
        mat *= norm
        generators.append(mat)
        
    return generators

def compute_and_print_d_values(N):
    """
    Computes and prints the unique non-zero d_ijk values for SU(N).
    """
    if N < 2:
        print(f"For SU({N}), the concept of d_ijk is trivial or undefined.")
        return 0

    lambdas = generate_gell_mann(N)
    num_gens = len(lambdas)
    d_values = set()
    
    print(f"Calculating d_ijk values for SU({N})...")

    # Iterate over all combinations of indices i, j, k with i <= j <= k
    # This is sufficient due to the total symmetry of d_ijk.
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                la, lb, lc = lambdas[i], lambdas[j], lambdas[k]
                
                # The formula is d_ijk = (1/4) * Tr({lambda_i, lambda_j} * lambda_k)
                # {A,B} = AB + BA
                d_val = 0.25 * np.trace((la @ lb + lb @ la) @ lc)
                
                # The trace of a product of Hermitian matrices is real.
                # We round to handle floating-point inaccuracies.
                d_rounded = np.round(d_val.real, 8)
                
                # Add non-zero values to the set.
                if abs(d_rounded) > 1e-7:
                    d_values.add(d_rounded)

    # Sort the values for a consistent and readable output
    sorted_values = sorted(list(d_values))
    
    print("-" * 30)
    print(f"For SU({N}), there are {len(sorted_values)} different non-zero d_ijk values.")
    
    if sorted_values:
        print("The unique values are:")
        # The prompt asks to "output each number in the final equation"
        # We interpret this as printing the distinct values found.
        for val in sorted_values:
            print(val)
    print("-" * 30)
            
    return len(sorted_values)

def main():
    """
    Main function to parse command-line arguments and run the calculation.
    """
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} N")
        print("       where N is an integer >= 2 for the group SU(N).")
        print("\nExample: python {sys.argv[0]} 3")
        return
        
    try:
        N = int(sys.argv[1])
        if N < 2:
            print("Error: N must be an integer greater than or equal to 2.")
            return
    except ValueError:
        print("Error: Invalid input. N must be an integer.")
        return
        
    compute_and_print_d_values(N)

if __name__ == "__main__":
    main()