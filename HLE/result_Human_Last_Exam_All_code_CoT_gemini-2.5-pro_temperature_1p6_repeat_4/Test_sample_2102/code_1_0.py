import numpy as np
from scipy.special import comb
from math import factorial

def solve():
    """
    Solves the problem by finding the smallest n where f(n) > 10 and then
    calculating n * ||W_n||_inf.
    """
    # Set a reasonable upper bound for n to search within.
    max_n_limit = 25 
    
    # Memoization tables for coefficients to speed up calculation.
    # a_k are coefficients for 2/pi * K(x), b_k are for e^x.
    a_coeffs = {}
    fact_cache = {i: factorial(i) for i in range(max_n_limit + 1)}

    def get_a_coeff(k):
        """Calculates and caches the coefficient a_k."""
        if k not in a_coeffs:
            # a_k = ((scipy.special.comb(2k, k) / 4**k))**2
            val = comb(2 * k, k, exact=True) / (4**k)
            a_coeffs[k] = val**2
        return a_coeffs[k]

    # Pre-populate a_k coefficients for efficiency.
    for k in range(max_n_limit + 1):
        get_a_coeff(k)

    # Loop through n to find the one satisfying the condition f(n) > 10.
    for n in range(1, max_n_limit + 1):
        # Calculate coefficients c_0, ..., c_n for the Taylor polynomial P_n(x).
        c_coeffs = []
        for i in range(n + 1):
            ci = sum(get_a_coeff(k) * (1 / fact_cache[i - k]) for k in range(i + 1))
            c_coeffs.append(ci)

        # The eigenvalues of the companion matrix are the roots of the polynomial.
        # numpy.roots expects coefficients from the highest power to the lowest.
        poly_coeffs_for_numpy = np.array(c_coeffs[::-1], dtype=float)
        
        eigenvalues = np.roots(poly_coeffs_for_numpy)

        # Calculate f(n) = sum of the absolute cubes of the eigenvalues.
        f_n = np.sum(np.abs(eigenvalues)**3)

        # Check if we have found the smallest n.
        if f_n > 10:
            # Determine the infinity norm of the Weyr matrix W_n.
            # This depends on the multiplicities of the eigenvalues.
            
            # Group close eigenvalues to find multiplicities due to floating point representation.
            # We round to a certain precision to count unique roots.
            rounded_eigenvalues = np.round(eigenvalues, decimals=8)
            unique_rounded_eigs, counts = np.unique(rounded_eigenvalues, return_counts=True)
            
            # Map rounded eigenvalues back to their original complex values to preserve precision.
            unique_original_eigs = []
            for ue in unique_rounded_eigs:
                # Find the first original eigenvalue that rounds to this unique value.
                idx = np.where(np.round(eigenvalues, decimals=8) == ue)[0][0]
                unique_original_eigs.append(eigenvalues[idx])

            # Calculate the norm of each Weyr block and find the maximum.
            max_norm = 0.0
            for i in range(len(unique_original_eigs)):
                eig_val = unique_original_eigs[i]
                count = counts[i]
                
                # For multiplicity > 1, the Weyr block norm is 1 + |lambda|.
                # For multiplicity 1, it's just |lambda|.
                norm = np.abs(eig_val)
                if count > 1:
                    norm += 1
                
                if norm > max_norm:
                    max_norm = norm
            
            # We've found n and calculated the norm, now compute the final result.
            result = n * max_norm
            
            # Print the final equation as requested.
            print(f"{n} * {max_norm} = {result}")
            
            # Return the final numerical answer in the required format
            global final_answer
            final_answer = result
            return

# Run the solver
final_answer = None
solve()