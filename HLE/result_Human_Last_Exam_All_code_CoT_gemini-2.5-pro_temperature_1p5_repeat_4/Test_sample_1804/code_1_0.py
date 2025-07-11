import sys

# Set the finite field F_q, where q is the prime F
F = 7

# --- Polynomial Arithmetic Functions over F_q ---

def poly_trim(p):
    """Removes leading zero coefficients from a polynomial represented as a list."""
    # The zero polynomial is represented as an empty list for convenience.
    if not p:
        return []
    p_copy = list(p)
    while len(p_copy) > 1 and p_copy[-1] == 0:
        p_copy.pop()
    if len(p_copy) == 1 and p_copy[0] == 0:
        return []
    return p_copy

def poly_degree(p):
    """Calculates the degree of a polynomial."""
    p_trimmed = poly_trim(p)
    return len(p_trimmed) - 1

def poly_add(p1, p2):
    """Adds two polynomials."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res[i] = (c1 + c2) % F
    return poly_trim(res)

def poly_sub(p1, p2):
    """Subtracts two polynomials."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res[i] = (c1 - c2 + F) % F
    return poly_trim(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials."""
    deg1 = poly_degree(p1)
    deg2 = poly_degree(p2)
    if deg1 == -1 or deg2 == -1:
        return []
    res_deg = deg1 + deg2
    res = [0] * (res_deg + 1)
    for i in range(len(p1)):
        for j in range(len(p2)):
            res[i+j] = (res[i+j] + p1[i] * p2[j]) % F
    return poly_trim(res)

def poly_divmod(N, D):
    """Performs polynomial division with remainder (N / D)."""
    degD = poly_degree(D)
    if degD == -1:
        raise ZeroDivisionError("Polynomial division by zero")

    N_rem = list(N)
    degN_rem = poly_degree(N_rem)
    
    if degN_rem < degD:
        return ([0], N_rem)

    q = [0] * (degN_rem - degD + 1)
    inv_d_lc = pow(D[-1], F - 2, F)
    
    while degN_rem >= degD:
        c = (N_rem[-1] * inv_d_lc) % F
        pos = degN_rem - degD
        q[pos] = c
        
        term_to_sub = [0] * pos + D
        term_to_sub = poly_mul(term_to_sub, [c])
        
        N_rem = poly_sub(N_rem, term_to_sub)
        degN_rem = poly_degree(N_rem)
        
    return poly_trim(q), N_rem

def poly_gcd(p1, p2):
    """Computes the greatest common divisor (GCD) of two polynomials."""
    a, b = p1, p2
    while poly_degree(b) != -1:
        _, r = poly_divmod(a, b)
        a, b = b, r
    # Make monic
    if poly_degree(a) > -1:
      lc = a[-1]
      inv_lc = pow(lc, F - 2, F)
      a = poly_mul(a, [inv_lc])
    return a

def poly_pow_mod(base, exp, mod_poly):
    """Computes (base^exp) mod mod_poly using modular exponentiation."""
    res = [1]
    
    _, base = poly_divmod(base, mod_poly)

    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base)
            _, res = poly_divmod(res, mod_poly)
        base = poly_mul(base, base)
        _, base = poly_divmod(base, mod_poly)
        exp //= 2
    return res

# --- Irreducibility Test ---

def is_irreducible(p):
    """
    Checks if polynomial p is irreducible over F_q using Ben-Or's algorithm.
    """
    n = poly_degree(p)
    if n <= 1:
        return True # By convention

    x = [0, 1]  # Polynomial x

    # A polynomial p(x) of degree n is irreducible iff
    # gcd(p(x), x^(q^k) - x) = 1 for k = 1, ..., floor(n/2)
    for k in range(1, n // 2 + 1):
        q_k = F**k
        
        # Calculate x^(q^k) mod p
        term = poly_pow_mod(x, q_k, p)
        
        # g = x^(q^k) - x
        g = poly_sub(term, x)
        
        # If gcd(p, g) is not 1 (degree > 0), then p is reducible
        d = poly_gcd(p, g)
        if poly_degree(d) > 0:
            return False
            
    return True

# --- Main Logic ---

def solve_problem():
    """
    Finds the set A and computes the final result.
    """
    A = []
    # Iterate through all possible values of 'a' in the field F
    for a_val in range(F):
        # Construct the polynomial p(x) = x^5 + a*x + 3
        # Represented as [c0, c1, c2, c3, c4, c5]
        # So, [3, a_val, 0, 0, 0, 1]
        poly = [3, a_val, 0, 0, 0, 1]
        
        if is_irreducible(poly):
            A.append(a_val)
            
    print(f"The set A of elements a for which x^5+ax+3 is irreducible is: {A}")
    
    if not A:
        print("Set A is empty. Cannot perform calculation.")
        return

    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    result = max_A ** min_A - len_A
    
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {len_A}")
    print("\nThe final calculation is:")
    print(f"{max_A}^{min_A} - {len_A} = {result}")
    
    # Required format for the final answer
    sys.stdout.write(f"<<<{result}>>>\n")

if __name__ == '__main__':
    solve_problem()