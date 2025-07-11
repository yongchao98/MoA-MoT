import math
import numpy as np

def solve_problem():
    """
    This function implements the plan to find the smallest n such that f(n) > 10
    and then calculates n * ||W_n||_inf.
    """
    # Cache for factorials and combinations to speed up calculations
    _fact_cache = {0: 1}
    def factorial(n):
        if n not in _fact_cache:
            _fact_cache[n] = n * factorial(n - 1)
        return _fact_cache[n]

    _comb_cache = {}
    def combinations(n, k):
        if k < 0 or k > n:
            return 0
        if (n, k) in _comb_cache:
            return _comb_cache[(n, k)]
        if k > n // 2:
            k = n - k
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        _comb_cache[(n, k)] = res
        return res

    # Caches for the series coefficients
    a_coeffs = []
    inv_factorials = []
    c_coeffs = []

    def get_a(k):
        nonlocal a_coeffs
        if k < len(a_coeffs):
            return a_coeffs[k]
        # Grow the cache if needed
        for i in range(len(a_coeffs), k + 1):
            val = (combinations(2 * i, i) / (4**i))**2
            a_coeffs.append(val)
        return a_coeffs[k]

    def get_inv_fact(k):
        nonlocal inv_factorials
        if k < len(inv_factorials):
            return inv_factorials[k]
        for i in range(len(inv_factorials), k + 1):
            inv_factorials.append(1.0 / factorial(i))
        return inv_factorials[k]

    # Main loop to find n
    n = 1
    while True:
        n += 1
        
        # Ensure c_coeffs cache is large enough for polynomial of degree n
        if len(c_coeffs) < n + 1:
            for m in range(len(c_coeffs), n + 1):
                c_m = 0.0
                for k in range(m + 1):
                    c_m += get_a(k) * get_inv_fact(m - k)
                c_coeffs.append(c_m)

        # Get coefficients for P_n(x) = c_0 + c_1*x + ... + c_n*x^n
        current_poly_coeffs = c_coeffs[:n+1]
        
        # numpy.roots expects coefficients in order [c_n, c_{n-1}, ..., c_0]
        poly_for_roots = np.array(current_poly_coeffs)[::-1].astype(float)
        
        # Find the roots (which are the eigenvalues of W_n)
        roots = np.roots(poly_for_roots)
        
        # Calculate f(n) = sum of absolute cubes of eigenvalues
        f_n = np.sum(np.abs(roots)**3)
        
        if f_n > 10:
            target_n = n
            
            # Check for multiple roots to calculate the infinity norm of W_n
            is_multiple = [False] * target_n
            # Use a tolerance for floating-point comparisons
            tolerance = 1e-8
            for i in range(target_n):
                for j in range(i + 1, target_n):
                    if np.isclose(roots[i], roots[j], atol=tolerance):
                        is_multiple[i] = True
                        is_multiple[j] = True
            
            # Calculate the norm contribution for each root
            # norm is |lambda| for simple roots, |lambda|+1 for multiple roots
            root_norms = []
            for i in range(target_n):
                norm = np.abs(roots[i])
                if is_multiple[i]:
                    norm += 1
                root_norms.append(norm)
            
            # ||W_n||_inf is the maximum of these individual root norms
            W_n_inf_norm = np.max(root_norms)
            
            # Calculate the final result
            final_answer = target_n * W_n_inf_norm
            
            # Print the components of the final equation as requested
            print(f"{target_n} * {W_n_inf_norm} = {final_answer}")
            
            # Print the final answer in the specified format
            print(f"<<<{final_answer}>>>")
            return

# Run the solver
solve_problem()