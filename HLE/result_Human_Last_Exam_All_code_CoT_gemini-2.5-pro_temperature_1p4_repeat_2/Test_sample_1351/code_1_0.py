import math
from fractions import Fraction

def get_divisors(n):
    """Returns all divisors of a positive integer n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return list(divs)

_mobius_cache = {}
def mobius(n):
    """Computes the Mobius function mu(n)."""
    if n in _mobius_cache:
        return _mobius_cache[n]
    if n == 1:
        return 1
    
    factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            count = 0
            while temp_n % d == 0:
                count += 1
                temp_n //= d
            if count > 1:
                _mobius_cache[n] = 0
                return 0
            factors.append(d)
        d += 1
    if temp_n > 1:
        factors.append(temp_n)
    
    _mobius_cache[n] = (-1)**len(factors)
    return _mobius_cache[n]

def N_k(k, q):
    """Computes the number of monic irreducible polynomials of degree k over F_q."""
    if k == 0: return 0
    divs = get_divisors(k)
    s = sum(mobius(j) * (q**(k//j)) for j in divs)
    return s // k

def product_term(start, end, func):
    """Computes the product of func(i) for i from start to end-1."""
    res = 1
    for i in range(start, end):
        res *= func(i)
    return res

def q_binom(n, k, q):
    """Computes the q-binomial coefficient."""
    if k < 0 or k > n:
        return 0
    num = product_term(0, k, lambda i: q**(n - i) - 1)
    den = product_term(1, k + 1, lambda i: q**i - 1)
    return num // den

def GL_size(k, q):
    """Computes the size of the general linear group GL_k(q)."""
    if k == 0: return 1
    return product_term(0, k, lambda i: q**k - q**i)

def I_k(k, q):
    """Computes the number of elements in GL_k(q) with an irreducible characteristic polynomial."""
    if k == 0: return 0
    num_irred_poly = N_k(k, q)
    glk_size = GL_size(k, q)
    return (num_irred_poly * glk_size) // (q**k - 1)

def solve():
    """Main calculation function."""
    d = 5
    q = 4
    e1 = 3
    e2 = 2

    # Compute the components of the formula
    q_binom_val = q_binom(d, e1, q)
    I_e1_val = I_k(e1, q)
    I_e2_val = I_k(e2, q)
    q_power_e1e2 = q**(e1 * e2)

    # Number of irreducible duos using the formula:
    # N_irred = q_binom(d, e1, q) * q**(e1*e2) * I_k(e1, q) * I_k(e2, q) * (q**(e1*e2) - 1)**2
    N_irred_duos = q_binom_val * q_power_e1e2 * I_e1_val * I_e2_val * (q_power_e1e2 - 1)**2

    # Total number of pairs in G x G
    G_size = GL_size(d, q)
    total_pairs = G_size**2

    # Proportion
    proportion = Fraction(N_irred_duos, total_pairs)

    print(f"(a) No\n(b) {{(2), (3)}}")
    print("(c) The proportion is calculated as N / |G|^2, where:")
    print(f"N (Number of irreducible (3,2)-stingray duos) = {N_irred_duos}")
    print(f"|G|^2 (Total pairs in GL_{d}({q}) x GL_{d}({q})) = {total_pairs}")
    print(f"Proportion = {float(proportion):.5e}")

solve()