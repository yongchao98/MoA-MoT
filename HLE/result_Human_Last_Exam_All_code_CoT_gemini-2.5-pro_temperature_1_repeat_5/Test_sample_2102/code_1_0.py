import math
import numpy as np
from math import lgamma

def solve():
    """
    This function implements the plan to find the final answer.
    """

    def get_taylor_coeffs(n_max):
        """
        Computes Taylor coefficients c_0, ..., c_{n_max} for (2/pi)K(x)e^x.
        Uses lgamma for numerical stability.
        Returns coeffs in order c_{n_max}, c_{n-1}, ..., c_0 for numpy.roots.
        """
        # Coefficients for (2/pi)K(x) = sum(a_k * x^k)
        # a_k = [binom(2k, k) / 4^k]^2
        a = []
        for k in range(n_max + 1):
            if k == 0:
                a.append(1.0)
            else:
                # Use lgamma for log(k!) to avoid overflow with large k
                # log(comb(n, k)) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
                log_comb_2k_k = lgamma(2 * k + 1) - 2 * lgamma(k + 1)
                log_a_k = 2 * (log_comb_2k_k - k * math.log(4))
                a.append(math.exp(log_a_k))

        # Coefficients for e^x = sum(b_k * x^k), b_k = 1/k!
        b_inv_fact = [1.0 / math.factorial(k) for k in range(n_max + 1)]

        # Coefficients for the product series (Cauchy Product)
        # c_i = sum_{k=0 to i} a_k * b_{i-k}
        c = []
        for i in range(n_max + 1):
            c_i = sum(a[k] * b_inv_fact[i - k] for k in range(i + 1))
            c.append(c_i)

        # numpy.roots expects coefficients of decreasing powers: [c_n, c_{n-1}, ..., c_0]
        return c[::-1]

    # Find the smallest n where f(n) > 10
    n = 1
    final_eigenvalues = None
    n_final = -1

    while True:
        n += 1
        # The Taylor polynomial of degree n has n+1 coefficients (c_0 to c_n)
        coeffs = get_taylor_coeffs(n)
        
        # The eigenvalues are the roots of the polynomial P_n(x)
        eigenvalues = np.roots(coeffs)
        
        # f(n) is the sum of the absolute cubes of the eigenvalues
        f_n = np.sum(np.abs(eigenvalues)**3)
        
        if f_n > 10:
            n_final = n
            final_eigenvalues = eigenvalues
            break
        
        if n > 100: # Safety break in case the condition is not met
            print("Error: Could not find n in a reasonable range.")
            return

    # Calculate n * ||W_n||_inf for the found n
    if n_final != -1:
        # Check for multiplicity of eigenvalues to determine the Jordan/Weyr structure.
        # We round to a high precision to check for numerically equal roots.
        precision = 8
        rounded_eigs = np.round(final_eigenvalues, precision)
        unique_eigs, counts = np.unique(rounded_eigs, return_counts=True)
        
        max_norm = 0.0
        # The infinity norm of the full Weyr matrix is the max of its block norms.
        for val, m in zip(unique_eigs, counts):
            # Infinity norm of a Weyr block for eigenvalue 'val' of size 'm'
            # is |val| if m=1, and |val|+1 if m > 1.
            block_norm = np.abs(val) + (1 if m > 1 else 0)
            if block_norm > max_norm:
                max_norm = block_norm
                
        W_n_inf_norm = max_norm
        
        result = n_final * W_n_inf_norm
        
        # As requested, output each number in the final equation.
        print(f"{n_final} * {W_n_inf_norm} = {result}")

solve()