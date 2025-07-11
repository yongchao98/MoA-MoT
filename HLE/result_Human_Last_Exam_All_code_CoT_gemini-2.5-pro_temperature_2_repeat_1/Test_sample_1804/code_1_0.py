import math

def solve():
    """
    Finds the set A and computes the final expression.
    """
    F = 7

    # --- Helper functions for polynomial arithmetic in F_7 ---
    # Polynomials are represented as lists of coefficients, e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2

    def poly_trim(p):
        """Removes leading zero coefficients."""
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p

    def poly_add(p1, p2):
        """Adds two polynomials in F_7."""
        n = max(len(p1), len(p2))
        p1_padded = p1 + [0] * (n - len(p1))
        p2_padded = p2 + [0] * (n - len(p2))
        res = [(c1 + c2) % F for c1, c2 in zip(p1_padded, p2_padded)]
        return poly_trim(res)

    def poly_mul(p1, p2):
        """Multiplies two polynomials in F_7."""
        n1, n2 = len(p1), len(p2)
        res = [0] * (n1 + n2 - 1)
        for i in range(n1):
            for j in range(n2):
                res[i + j] = (res[i + j] + p1[i] * p2[j]) % F
        return poly_trim(res)
    
    def poly_divmod(dividend, divisor):
        """Performs polynomial long division in F_7. Returns (quotient, remainder)."""
        mod_inv = {1: 1, 2: 4, 3: 5, 4: 2, 5: 3, 6: 6}
        rem = list(dividend)
        divisor = poly_trim(divisor)
        
        if divisor == [0]:
            raise ZeroDivisionError("Polynomial division by zero")
            
        deg_d = len(divisor) - 1
        
        if len(rem) < len(divisor):
            return [0], rem

        quot = [0] * (len(rem) - deg_d)
        lead_d_inv = mod_inv[divisor[-1]]

        for i in range(len(rem) - 1, deg_d - 1, -1):
            if i < len(rem):
                coeff = (rem[i] * lead_d_inv) % F
                deg_q = i - deg_d
                quot[deg_q] = coeff
                for j in range(deg_d + 1):
                    term_index = deg_q + j
                    sub_val = (coeff * divisor[j]) % F
                    rem[term_index] = (rem[term_index] - sub_val + F) % F
        
        return poly_trim(quot), poly_trim(rem)

    def poly_gcd(p1, p2):
        """Computes the GCD of two polynomials in F_7."""
        a, b = p1, p2
        while b != [0]:
            _, r = poly_divmod(a, b)
            a, b = b, r
        return a
        
    def poly_pow_mod(base, exp, modulus):
        """Computes (base^exp) % modulus for polynomials."""
        res = [1]
        current_power = list(base)
        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, current_power)
                _, res = poly_divmod(res, modulus)
            current_power = poly_mul(current_power, current_power)
            _, current_power = poly_divmod(current_power, modulus)
            exp //= 2
        return res

    def is_irreducible_deg5(p):
        """
        Checks if a polynomial of degree 5 is irreducible over F_7 using GCD tests.
        """
        # Test 1: Check for roots in F_7. This is faster than a full GCD.
        has_root = False
        for k in range(F):
            val = (pow(k, 5, F) + p[1] * k + p[0]) % F
            if val == 0:
                has_root = True
                break
        if has_root:
            return False

        # Test 2: Check for irreducible quadratic factors via gcd(p, x^49 - x).
        # We compute h(x) = x^49 mod p(x), then find gcd(h(x) - x, p(x)).
        x_poly = [0, 1]
        h_x = poly_pow_mod(x_poly, F * F, p)
        h_x_minus_x = poly_add(h_x, [0, F - 1]) # Subtract x
        gcd_res = poly_gcd(h_x_minus_x, p)

        # The polynomial is irreducible if the gcd is a constant (degree 0).
        return len(gcd_res) == 1

    A = []
    for a in range(F):
        # p(x) = x^5 + ax + 3 is represented as [3, a, 0, 0, 0, 1]
        poly = [3, a, 0, 0, 0, 1]
        if is_irreducible_deg5(poly):
            A.append(a)

    print(f"The set A of coefficients 'a' for which the polynomial is irreducible is: {A}")
    
    if not A:
        print("The set A is empty. Cannot perform calculation.")
        return

    min_A = min(A)
    max_A = max(A)
    card_A = len(A)
    
    # Perform the final calculation using standard integer arithmetic
    result = int(math.pow(max_A, min_A) - card_A)

    print(f"The minimum element of A is: {min_A}")
    print(f"The maximum element of A is: {max_A}")
    print(f"The size of A is: {card_A}")
    
    print(f"\nThe final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {card_A} = {result}")

solve()
<<<3>>>