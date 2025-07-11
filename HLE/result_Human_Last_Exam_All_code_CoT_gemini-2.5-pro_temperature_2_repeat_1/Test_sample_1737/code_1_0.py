import numpy as np
import math

def get_sun_generators(N):
    """
    Constructs a standard basis of generators for the SU(N) Lie algebra.
    The generators are N x N traceless Hermitian matrices, normalized
    such that Tr(T^a T^b) = (1/2) * delta_ab.
    """
    if N <= 1:
        return []

    generators = []
    
    # Symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            S = np.zeros((N, N), dtype=float)
            S[j, k] = 1.0
            S[k, j] = 1.0
            # Normalization: Tr(S^2/c^2) = 1/2 -> (1/c^2) * Tr(diag(1,1,0..)) = 1/2 -> 2/c^2=1/2 -> c=2
            generators.append(S / 2.0)

    # Anti-symmetric off-diagonal generators
    for j in range(N):
        for k in range(j + 1, N):
            A = np.zeros((N, N), dtype=complex)
            A[j, k] = -1.0j
            A[k, j] = 1.0j
            # Normalization is the same as symmetric case
            generators.append(A / 2.0)
            
    # Diagonal generators
    for l_idx in range(1, N):
        diag_vals = np.zeros(N, dtype=float)
        diag_vals[:l_idx] = 1.0
        diag_vals[l_idx] = -l_idx
        
        # Normalization: Tr(D^2/c^2) = 1/2 -> (1/c^2) * sum(diag_vals^2) = 1/2
        # sum(diag_vals^2) = l_idx*1^2 + (-l_idx)^2 = l_idx + l_idx^2 = l_idx*(l_idx+1)
        norm_factor = 1.0 / np.sqrt(2 * l_idx * (l_idx + 1))
        D = np.diag(diag_vals * norm_factor)
        generators.append(D)
        
    return generators

def solve():
    """
    Calculates the number of distinct non-zero values for the totally symmetric
    structure constants d_ijk of SU(N).
    """
    # Set the value for N here.
    N = 3

    print(f"Calculating for SU({N})...")

    if N <= 1:
        print("For SU(N) with N<=1, the Lie algebra is trivial.")
        print("Number of distinct non-zero d_ijk values: 0")
        return
    
    if N == 2:
        print("For SU(2), all d_ijk constants are zero.")
        print("Number of distinct non-zero d_ijk values: 0")
        return

    generators = get_sun_generators(N)
    num_generators = len(generators)
    
    # A small number for floating point comparisons
    epsilon = 1e-9
    
    # Use a set to store unique d_ijk values, rounded to avoid float precision issues.
    d_values = set()
    
    # Loop over all combinations with i <= j <= k to exploit symmetry
    for i in range(num_generators):
        for j in range(i, num_generators):
            # Pre-calculate the anti-commutator {T_i, T_j}
            ti = generators[i]
            tj = generators[j]
            anti_commutator = ti @ tj + tj @ ti
            
            for k in range(j, num_generators):
                tk = generators[k]
                
                # Calculate d_ijk = 2 * Tr({T_i, T_j} * T_k)
                d_ijk = 2 * np.trace(anti_commutator @ tk)
                
                # Result should be real; take real part to discard small imaginary noise
                d_val = d_ijk.real
                
                if abs(d_val) > epsilon:
                    # Round to handle floating point variations
                    rounded_val = round(d_val, 8)
                    d_values.add(rounded_val)

    # Output the result
    print(f"\nThe non-zero d_ijk for SU({N}) take the following {len(d_values)} distinct values:")
    if d_values:
        # Sort for consistent display
        sorted_values = sorted(list(d_values))
        print(sorted_values)
    else:
        print("None")
        
    print(f"\nFinal Answer: The total number of distinct non-zero values is {len(d_values)}")

if __name__ == '__main__':
    solve()