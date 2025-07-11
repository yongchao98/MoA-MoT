import numpy as np
from scipy.special import factorial, binom
from collections import Counter

def solve():
    """
    Solves the problem by finding the smallest n where f(n) > 10
    and then calculating n * ||W_n||_inf.
    """
    # Memoization caches for calculated coefficients to speed up computation
    memo_a = {}
    memo_c = {}

    def get_a(k):
        """
        Calculates the k-th coefficient of the Taylor series for 2/pi * K(x).
        a_k = (binom(2k, k) / 4^k)^2
        """
        if k in memo_a:
            return memo_a[k]
        if k == 0:
            val = 1.0
        else:
            # Using scipy.special.binom which handles large numbers and returns a float
            val = (binom(2 * k, k) / (4**k))**2
        memo_a[k] = val
        return val

    def get_c(n):
        """
        Calculates the n-th coefficient of the Taylor series for g(x) = 2/pi * K(x) * e^x
        using the Cauchy product of the series for 2/pi*K(x) and e^x.
        """
        if n in memo_c:
            return memo_c[n]
        # c_n = sum_{k=0 to n} a_k * b_{n-k} where b_j = 1/j!
        c_n = sum(get_a(k) / factorial(n - k, exact=False) for k in range(n + 1))
        memo_c[n] = c_n
        return c_n

    target_n = -1
    eigenvalues = None
    f_n_value = 0
    max_n = 50  # A reasonable limit to prevent an infinite loop

    for n in range(1, max_n + 1):
        # The coefficients for numpy.roots must be in order of descending power
        coeffs = [get_c(k) for k in range(n, -1, -1)]
        
        # The eigenvalues of S_n are the roots of the Taylor polynomial
        roots = np.roots(coeffs)
        
        # f(n) is the sum of the absolute cubes of the eigenvalues
        f_n = np.sum(np.abs(roots)**3)
        
        if f_n > 10:
            target_n = n
            eigenvalues = roots
            f_n_value = f_n
            break

    if target_n == -1:
        print(f"No solution found for n up to {max_n}.")
        return

    # To check for repeated roots, we round them to a certain precision
    # and count the occurrences of each unique rounded root.
    precision = 8
    root_map = {}
    for r in eigenvalues:
        # A tuple of the rounded real and imaginary parts can serve as a hashable key
        rounded_r_key = (round(r.real, precision), round(r.imag, precision))
        if rounded_r_key not in root_map:
            root_map[rounded_r_key] = {'original': r, 'count': 0}
        root_map[rounded_r_key]['count'] += 1

    # Calculate the infinity norm of W_n
    # The norm of a Jordan block for eigenvalue λ with multiplicity m > 1 is |λ|+1.
    # If m=1, the norm is |λ|.
    max_norm_val = 0.0
    for data in root_map.values():
        original_root = data['original']
        multiplicity = data['count']
        
        current_norm = np.abs(original_root)
        if multiplicity > 1:
            current_norm += 1.0
        
        if current_norm > max_norm_val:
            max_norm_val = current_norm
    
    infinity_norm_Wn = max_norm_val
    result = target_n * infinity_norm_Wn

    print(f"The smallest integer n where f(n) > 10 is found to be n = {target_n}.")
    print(f"For n = {target_n}, the value f(n) is {f_n_value:.4f}.")
    print(f"The eigenvalues for n = {target_n} appear to be distinct.")
    print(f"Thus, the infinity norm ||W_{target_n}||_inf is the maximum absolute value of the eigenvalues.")
    print(f"The value of the infinity norm ||W_{target_n}||_inf is {infinity_norm_Wn:.4f}.")
    print(f"The final result is n * ||W_n||_inf = {target_n} * {infinity_norm_Wn:.4f} = {result:.4f}.")
    print(f"<<<{result}>>>")

solve()