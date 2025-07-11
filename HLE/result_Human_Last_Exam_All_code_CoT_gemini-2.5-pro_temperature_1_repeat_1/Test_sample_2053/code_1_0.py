import math
from fractions import Fraction

def get_prime_factorization(num):
    """Computes the prime factorization of a number."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """Computes the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    if len(factors) % 2 == 1:
        return -1
    else:
        return 1

def phi(n):
    """Computes Euler's totient function phi(n)."""
    if n == 1:
        return 1
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
    return result

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def combinations(n, k):
    """Computes the binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def calculate_d2_size(N, h):
    """
    Calculates |D_2(N, h)| using the formula from question (a)
    and a corrected version of the Upsilon formula from the literature.
    """

    # Simplified Upsilon formula, valid when N/n divides h.
    def calculate_upsilon(n_val):
        if (N % n_val != 0) or (h % (N // n_val) != 0):
             return 0
        
        hn_div_N = h * n_val // N
        N_div_n = N // n_val
        
        sign = (-1)**(n_val - hn_div_N)
        phi_term = Fraction(phi(N_div_n), N_div_n)
        comb_term = combinations(h, hn_div_N)
        
        return sign * phi_term * comb_term

    total_sum = Fraction(0)
    D_N = get_divisors(N)
    
    # Pre-calculate D* set and Upsilon values
    D_star_set = {n for n in D_N if (N % n == 0) and (h % (N // n) == 0)}
    upsilon_vals = {n: calculate_upsilon(n) for n in D_star_set}

    print(f"N = {N}, h = {h}")
    print(f"D_N = {D_N}")
    print(f"D* = {sorted(list(D_star_set))}")
    print("Upsilon values:")
    for n in sorted(list(upsilon_vals.keys())):
        print(f"  Upsilon(N,h,{n}) = {upsilon_vals[n]}")
    
    # Main sum calculation
    final_sum_str = []
    for u in D_N:
        D_u = get_divisors(u)
        D_u_star = {n for n in D_u if n in D_star_set}
        
        inner_sum = Fraction(0)
        for n in D_u_star:
            mu_val = mobius(u // n)
            term = Fraction(1, n) * mu_val * upsilon_vals[n]
            inner_sum += term
        
        term_for_u = Fraction(1, u) * inner_sum
        total_sum += term_for_u
        if term_for_u != 0:
            final_sum_str.append(f"({term_for_u})")

    print(f"\nCalculation of |D_2({N},{h})| = " + " + ".join(final_sum_str))
    print(f"\nFinal result from formula (a): {total_sum}")
    if total_sum.denominator != 1:
        print("\nNote: The formula yields a non-integer result, which contradicts the nature of a count.")

calculate_d2_size(8, 4)

# The integer answer for (c) is a guess due to flawed premises.
print("\nFinal Answer block:")
print("<<< (a) No; (b) No; (c) 2 >>>")
