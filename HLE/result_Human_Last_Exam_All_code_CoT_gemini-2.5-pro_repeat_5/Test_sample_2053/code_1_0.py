import math
from fractions import Fraction

# --- Memoization dictionaries ---
memo_phi = {}
memo_mu = {}
memo_divs = {}
memo_Sigma = {}
memo_Upsilon = {}

def phi(n):
    """Computes Euler's totient function."""
    if n in memo_phi:
        return memo_phi[n]
    result = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            result -= result // p
        p += 1
    if temp_n > 1:
        result -= result // temp_n
    memo_phi[n] = result
    return result

def mu(n):
    """Computes the Möbius function."""
    if n in memo_mu:
        return memo_mu[n]
    if n == 1:
        memo_mu[n] = 1
        return 1
    
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            count = 0
            while temp_n % d == 0:
                count += 1
                temp_n //= d
            if count > 1:
                memo_mu[n] = 0
                return 0
            factors[d] = count
        d += 1
    if temp_n > 1:
        factors[temp_n] = 1
        
    if any(count > 1 for count in factors.values()):
        memo_mu[n] = 0
        return 0
    
    result = (-1)**len(factors)
    memo_mu[n] = result
    return result

def divisors(n):
    """Finds all positive divisors of an integer."""
    if n in memo_divs:
        return memo_divs[n]
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    result = sorted(list(divs))
    memo_divs[n] = result
    return result

def Sigma(j, n):
    """Computes the Sigma_j^(n) term."""
    if (j, n) in memo_Sigma:
        return memo_Sigma[(j, n)]
    
    if not (0 <= j <= n):
        raise ValueError("j must be between 0 and n")
    if j == n:
        return Fraction(0)
    if j == 0:
        return Fraction(math.factorial(n - 1) - 1)

    # Sum part for j in [1, n-1]
    sum_val = Fraction(0)
    for k in range(n - j):
        denominator = math.factorial(k) * (j + k) * (n - j - k)
        if denominator == 0: continue
        term = Fraction((-1)**k, denominator)
        sum_val += term
    
    first_part = Fraction(math.factorial(n), math.factorial(j - 1)) * sum_val
    second_part = (-1)**(n - j) * math.comb(n - 1, j - 1)
    
    result = first_part + Fraction(second_part) - Fraction(1)
    memo_Sigma[(j, n)] = result
    return result

def power(base, exp):
    """Custom power function to handle 0^0=1 case."""
    if base == 0 and exp == 0:
        return 1
    return base**exp

def Upsilon(N, h, n):
    """Computes the Upsilon_{N,h,n} term."""
    if (N, h, n) in memo_Upsilon:
        return memo_Upsilon[(N, h, n)]

    # From the definition of D*, hn/N is an integer.
    hn_div_N = (h * n) // N
    
    # Pre-calculate the term B, as it doesn't depend on m
    n_prime = n - hn_div_N
    n_div_N = Fraction(n, N)
    n_minus_hn_div_N = n - hn_div_N
    
    powB1_base = n_div_N - 1
    powB1 = power(powB1_base, n_prime)
    powB2 = power(Fraction(-1), n_minus_hn_div_N)
    
    term_B_inner = n_div_N * (powB1 - powB2) + powB2
    term_B = Fraction(hn_div_N) * term_B_inner

    # The sum part (A)
    sum_val = Fraction(0)
    m_start = hn_div_N # Since hn/N is an integer
    for m in range(m_start, n):
        # Calculate term A_m
        phi_val = phi(N // n)
        
        pow1 = Fraction(N**(n - m - 1), n**(n - m - 1))
        
        base2 = Fraction(m, N)
        exp2 = m - hn_div_N
        pow2 = power(base2, exp2)

        factor_1n = Fraction(1, n) - 1
        
        sigma_diff = Sigma(m, n) - Sigma(m + 1, n)
        
        term_A_m = phi_val * pow1 * pow2 * factor_1n * sigma_diff
        
        # The summand is (A_m + B)
        sum_val += (term_A_m + term_B)

    memo_Upsilon[(N, h, n)] = sum_val
    return sum_val

def calculate_D2(N, h):
    """Calculates |D_2(N, h)|."""
    D_N = divisors(N)
    
    # D* = {n in D_N; N/n | h}
    D_star = {n for n in D_N if h % (N // n) == 0}

    total_sum = Fraction(0)
    
    # The formula can be simplified by collecting coefficients of Upsilon
    # Coeff(n) = sum_{u | N, n | u} (1/u) * (1/n) * mu(u/n)
    upsilon_coeffs = {}
    
    for u in D_N:
        D_u = divisors(u)
        D_u_star = [n for n in D_u if n in D_star]
        
        for n in D_u_star:
            mu_val = mu(u // n)
            if mu_val == 0:
                continue
            
            coeff = Fraction(1, u) * Fraction(1, n) * mu_val
            if n not in upsilon_coeffs:
                upsilon_coeffs[n] = Fraction(0)
            upsilon_coeffs[n] += coeff

    print(f"Calculating |D_2({N}, {h})|:")
    print("The formula is sum over n in D_star of Coeff(n) * Upsilon(N,h,n)")
    print("---------------------------------------------------------------")
    
    for n in sorted(upsilon_coeffs.keys()):
        coeff = upsilon_coeffs[n]
        upsilon_val = Upsilon(N, h, n)
        term = coeff * upsilon_val
        total_sum += term
        print(f"n = {n}: Coeff = {coeff}, Upsilon = {upsilon_val}")
        print(f"   Term contribution: {coeff} * {upsilon_val} = {term}")

    print("---------------------------------------------------------------")
    print(f"Final Sum (total |D_2({N}, {h})|) = {total_sum}")
    return total_sum

if __name__ == '__main__':
    N_val = 8
    h_val = 4
    
    result = calculate_D2(N_val, h_val)
    
    print("\nFinal Answer:")
    print(f"(a) Yes; (b) Yes; (c) {result.numerator}")
    print(f"<<<{result.numerator}>>>")
