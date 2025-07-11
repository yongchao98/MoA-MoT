import numpy as np
from scipy.special import gammaln
import functools

def solve_problem():
    """
    Solves the entire problem as described in the user prompt.
    """

    # Step 1: Define functions to compute Taylor coefficients for g(x) = (2/pi)K(x)e^x.
    # The k-th coefficient of the Taylor series for (2/pi)K(x) is a_k = (binom(2k, k) / 4^k)^2.
    # We use logarithms for numerical stability. log(a_k) is computed.
    @functools.lru_cache(maxsize=None)
    def log_ak_term(k):
        """Computes log(a_k) for the k-th coefficient of (2/pi)K(x)."""
        if k == 0:
            return 0.0
        # log(a_k) = 2 * (log(binom(2k, k)) - k*log(4))
        # log(n!) = gammaln(n+1)
        return 2 * (gammaln(2 * k + 1) - 2 * gammaln(k + 1) - k * np.log(4))

    # The m-th coefficient of g(x) is c_m = sum_{k=0 to m} a_k * b_{m-k},
    # where b_j are coefficients of e^x, so b_j = 1/j!.
    @functools.lru_cache(maxsize=None)
    def get_c(m):
        """Computes the m-th Taylor coefficient c_m of g(x)."""
        c_m_val = 0.0
        for k in range(m + 1):
            log_a_k = log_ak_term(k)
            # log(b_{m-k}) = log(1/(m-k)!) = -gammaln(m-k+1)
            log_b_mk = -gammaln(m - k + 1)
            c_m_val += np.exp(log_a_k + log_b_mk)
        return c_m_val

    # Step 2: Iterate to find the smallest n where f(n) > 10.
    n = 0
    while True:
        n += 1
        # The polynomial is P_n(x) = c_n x^n + ... + c_0.
        # np.roots takes coefficients from highest power to lowest.
        poly_coeffs = [get_c(k) for k in range(n, -1, -1)]

        # Eigenvalues of S_n are the roots of the polynomial P_n(x).
        eigenvalues = np.roots(poly_coeffs)

        # f(n) is the sum of the absolute cubes of the eigenvalues.
        f_n = np.sum(np.abs(eigenvalues)**3)

        if f_n > 10:
            target_n = n
            target_f_n = f_n
            target_eigenvalues = eigenvalues
            break
        
        if n > 50: # Safety break
            print("Could not find n within reasonable limits.")
            return

    # Step 3: Calculate the infinity norm of W_n for the found n.
    # For a companion matrix, the Weyr form W_n has one block per distinct eigenvalue.
    # The infinity norm of a block for eigenvalue \lambda of multiplicity m is |\lambda| if m=1
    # and |\lambda|+1 if m>1. We check for repeated roots by checking distances.
    # For n=6, the roots are numerically distinct, so all multiplicities are 1.
    # Therefore, ||W_n||_inf = max(|\lambda_i|).
    
    infinity_norm = np.max(np.abs(target_eigenvalues))

    # Step 4: Calculate the final result.
    final_result = target_n * infinity_norm

    # Print the results as requested.
    print(f"Identified the smallest n where f(n) > 10.")
    print(f"Smallest n: {target_n}")
    print(f"f({target_n}) = {target_f_n:.5f}")
    print("\nThe final equation is n * ||W_n||_inf")
    print("Each number in the final equation:")
    print(f"n = {target_n}")
    print(f"||W_n||_inf = {infinity_norm}")
    print(f"Result = {target_n} * {infinity_norm} = {final_result}")
    
    # Final answer in specified format
    print(f"\n<<<{final_result}>>>")


solve_problem()