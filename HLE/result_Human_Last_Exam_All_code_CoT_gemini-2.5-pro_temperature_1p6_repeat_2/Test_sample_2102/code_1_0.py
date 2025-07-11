import numpy as np
import math
from scipy.special import comb

def solve_and_print():
    """
    This function implements the plan to solve the problem.
    It finds the smallest integer n for which f(n) > 10, and then computes
    and prints the final expression n * ||W_n||_inf.
    """

    # Use memoization (caching) for coefficients to avoid recomputing them in each loop.
    c_coeffs_cache = {}
    a_coeffs_cache = {}

    def get_c_k(k):
        """Calculates (or retrieves from cache) the coefficient c_k for the series of 2/pi * K(x)."""
        if k in c_coeffs_cache:
            return c_coeffs_cache[k]
        
        # The formula for c_k is (binom(2k, k) / 4^k)^2
        # Use floating-point division as precision is sufficient.
        val = (comb(2 * k, k, exact=False) / (4**k))**2
        c_coeffs_cache[k] = val
        return val

    def get_a_n(n):
        """Calculates (or retrieves from cache) the coefficient a_n for the series of g(x)."""
        if n in a_coeffs_cache:
            return a_coeffs_cache[n]
        
        # The formula for a_n is the Cauchy product: sum_{k=0 to n} c_k / (n-k)!
        val = 0.0
        for k in range(n + 1):
            val += get_c_k(k) / math.factorial(n - k)
        
        a_coeffs_cache[n] = val
        return val

    # Initialize base cases for the caches
    c_coeffs_cache[0] = 1.0
    a_coeffs_cache[0] = 1.0

    n = 1
    while True:
        # Construct the list of polynomial coefficients for P_n(x), from highest power (a_n) to lowest (a_0).
        # np.roots expects this order.
        poly_coeffs = [get_a_n(i) for i in range(n, -1, -1)]
        
        # The eigenvalues of W_n are the roots of the polynomial P_n(x).
        eigenvalues = np.roots(poly_coeffs)
        
        # f(n) is the sum of the absolute cubes of the eigenvalues.
        f_n = np.sum(np.abs(eigenvalues)**3)
        
        # Check if we have found the smallest n such that f(n) > 10.
        if f_n > 10:
            # For the found n, calculate the infinity norm of W_n.
            # Assuming simple roots (which is true for this n), the Weyr matrix is
            # effectively diagonal, so its norm is the maximum absolute eigenvalue.
            inf_norm_Wn = np.max(np.abs(eigenvalues))
            
            # Calculate the final result.
            result = n * inf_norm_Wn
            
            # Print the final result in the requested format, showing each number.
            print(f"{n} * {inf_norm_Wn} = {result}")
            
            # Exit the loop as the solution is found.
            break
        
        n += 1

# Run the solver.
solve_and_print()