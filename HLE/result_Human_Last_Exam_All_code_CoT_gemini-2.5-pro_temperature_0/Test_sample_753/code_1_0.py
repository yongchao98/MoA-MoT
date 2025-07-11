import math
from collections import defaultdict
from functools import lru_cache

# Helper for partitions and cycle structures
@lru_cache(maxsize=None)
def get_partitions(n):
    if n == 0:
        return [[]]
    parts = []
    for i in range(1, n + 1):
        for p in get_partitions(n - i):
            if not p or i <= p[0]:
                parts.append([i] + p)
    return parts

def to_counts(p, n):
    counts = defaultdict(int)
    for part in p:
        counts[part] += 1
    return tuple(counts.get(i, 0) for i in range(1, n + 1))

def z_val(p_counts, n):
    val = 1
    for i in range(1, n + 1):
        c_i = p_counts[i-1]
        if c_i > 0:
            val *= (i**c_i) * math.factorial(c_i)
    return val

# Polynomials represented as dicts from exponent_tuple to coefficient
def poly_add(p1, p2):
    res = p1.copy()
    for k, v in p2.items():
        res[k] = res.get(k, 0) + v
    return res

def poly_scale(p, s):
    return {k: v * s for k, v in p.items()}

def poly_mult(p1, p2, n_vars):
    res = defaultdict(int)
    if not p1: return p2
    if not p2: return p1
    for k1, v1 in p1.items():
        for k2, v2 in p2.items():
            new_k = tuple(i + j for i, j in zip(k1, k2))
            res[new_k] += v1 * v2
    return dict(res)

def poly_pow(p, exp, n_vars):
    res = {tuple([0]*n_vars): 1}
    if exp == 0: return res
    base = p.copy()
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mult(res, base, n_vars)
        base = poly_mult(base, base, n_vars)
        exp //= 2
    return res

# Cycle index calculations
@lru_cache(maxsize=None)
def phi(n):
    res = n
    p = 2
    temp_n = n
    while p * p <= temp_n:
        if temp_n % p == 0:
            while temp_n % p == 0:
                temp_n //= p
            res -= res // p
        p += 1
    if temp_n > 1:
        res -= res // temp_n
    return res

@lru_cache(maxsize=None)
def get_z_cyclic(i, n_vars):
    poly = defaultdict(int)
    for d in range(1, i + 1):
        if i % d == 0:
            exp = [0] * n_vars
            exp[d-1] = i // d
            poly[tuple(exp)] += phi(d)
    return poly_scale(dict(poly), 1/i)

@lru_cache(maxsize=None)
def get_z_symmetric(c, n_vars):
    poly = defaultdict(int)
    for p in get_partitions(c):
        p_counts = to_counts(p, c)
        term_poly = {tuple([0]*n_vars): 1/z_val(p_counts, c)}
        for i in range(1, c + 1):
            if p_counts[i-1] > 0:
                exp = [0] * n_vars
                exp[i-1] = p_counts[i-1]
                monomial = {tuple(exp): 1}
                term_poly = poly_mult(term_poly, monomial, n_vars)
        poly = poly_add(poly, term_poly)
    return dict(poly)

def plethysm(P, Q, n_vars):
    res = {tuple([0]*n_vars): 0}
    for p_exp, p_coeff in P.items():
        term_res = {tuple([0]*n_vars): p_coeff}
        for i in range(n_vars):
            if p_exp[i] > 0:
                Q_i = defaultdict(int)
                for q_exp, q_coeff in Q.items():
                    new_exp = [0] * n_vars
                    for j in range(n_vars):
                        if q_exp[j] > 0:
                            if (j + 1) * (i + 1) -1 < n_vars:
                                new_exp[(j + 1) * (i + 1) - 1] = q_exp[j]
                    Q_i[tuple(new_exp)] = q_coeff
                
                term_res = poly_mult(term_res, poly_pow(dict(Q_i), p_exp[i], n_vars), n_vars)
        res = poly_add(res, term_res)
    return res

def solve():
    m = 3
    n = sum(k * (m + 1 - k) for k in range(1, m + 1))
    
    lambda_parts = [3, 2, 2, 1, 1, 1]
    lambda_counts = to_counts(lambda_parts, n)
    
    all_mu_parts = get_partitions(n)
    
    total_sum = 0
    
    print(f"Calculating cardinality of |Theta^-1(lambda)| for n={n} and lambda={lambda_parts}")
    print("-" * 80)
    print(f"{'Partition mu':<25} | {'N(lambda, mu)':>15} | {'(z_mu)^2':>25} | {'Contribution':>30}")
    print("-" * 80)

    for mu_parts in all_mu_parts:
        mu_counts = to_counts(mu_parts, n)
        
        # Calculate N(lambda, mu)
        centralizer_z_poly = {tuple([0]*n): 1}
        for i in range(1, n + 1):
            c_i = mu_counts[i-1]
            if c_i > 0:
                z_C_i = get_z_cyclic(i, n)
                z_S_ci = get_z_symmetric(c_i, n)
                wreath_z = plethysm(z_S_ci, z_C_i, n)
                centralizer_z_poly = poly_mult(centralizer_z_poly, wreath_z, n)
        
        z_mu_val = z_val(mu_counts, n)
        
        # Number of elements of type lambda in C(mu)
        # is |C(mu)| * coefficient of p_lambda in Z(C(mu))
        # p_lambda corresponds to the monomial for lambda
        p_lambda_exp = lambda_counts
        coeff = centralizer_z_poly.get(p_lambda_exp, 0)
        
        N_lambda_mu = round(z_mu_val * coeff)

        if N_lambda_mu > 0:
            z_mu_sq = z_mu_val**2
            contribution = N_lambda_mu * z_mu_sq
            total_sum += contribution
            
            mu_str = str(mu_parts)
            print(f"{mu_str:<25} | {N_lambda_mu:>15} | {z_mu_sq:>25} | {contribution:>30}")

    print("-" * 80)
    print(f"Total cardinality: {total_sum}")
    
    result_str = str(total_sum)
    print(f"\nFirst 40 digits of the result: {result_str[:40]}")
    
solve()
<<<1313837398313981319461120320335335321600>>>