def solve_polynomial_problem():
    """
    Solves the problem by finding the set A of values 'a' for which
    x^5 + ax + 3 is irreducible over F_7, and then calculating
    max(A)^min(A) - |A|.
    """
    
    # Define the modulus for the finite field F_7
    m = 7

    # Helper functions for polynomial arithmetic over F_m
    # Polynomials are represented as lists of coefficients, e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2.

    def normalize(p):
        """Removes leading zero coefficients from a polynomial representation."""
        p_copy = list(p)
        while len(p_copy) > 1 and p_copy[-1] == 0:
            p_copy.pop()
        return p_copy

    def poly_add(p1, p2):
        """Adds two polynomials over F_m."""
        n = max(len(p1), len(p2))
        res = [0] * n
        for i in range(len(p1)):
            res[i] = p1[i]
        for i in range(len(p2)):
            res[i] = (res[i] + p2[i]) % m
        return normalize(res)

    def poly_sub(p1, p2):
        """Subtracts the second polynomial from the first over F_m."""
        n = max(len(p1), len(p2))
        res = [0] * n
        for i in range(len(p1)):
            res[i] = p1[i]
        for i in range(len(p2)):
            res[i] = (res[i] - p2[i]) % m
        return normalize(res)

    def poly_mul(p1, p2):
        """Multiplies two polynomials over F_m."""
        n = len(p1) + len(p2) - 1
        res = [0] * n
        for i in range(len(p1)):
            for j in range(len(p2)):
                res[i+j] = (res[i+j] + p1[i] * p2[j]) % m
        return normalize(res)

    def mod_inverse(n):
        """Calculates the modular inverse of n modulo m (for prime m)."""
        return pow(n, m - 2, m)

    def poly_divmod(p1, p2):
        """Performs polynomial long division over F_m."""
        p1 = normalize(list(p1))
        p2 = normalize(list(p2))

        if len(p2) == 1 and p2[0] == 0:
            raise ZeroDivisionError("Polynomial division by zero")

        if len(p1) < len(p2):
            return [0], p1
        
        den_lead_inv = mod_inverse(p2[-1])
        q = [0] * (len(p1) - len(p2) + 1)
        
        while len(p1) >= len(p2):
            deg_diff = len(p1) - len(p2)
            lead_coeff = (p1[-1] * den_lead_inv) % m
            q[deg_diff] = lead_coeff
            
            term = [0] * (deg_diff + 1)
            term[-1] = lead_coeff
            
            sub = poly_mul(p2, term)
            p1 = poly_sub(p1, sub)
        
        return normalize(q), p1

    def poly_gcd(p1, p2):
        """Calculates the GCD of two polynomials over F_m."""
        while p2 != [0]:
            _, r = poly_divmod(p1, p2)
            p1, p2 = p2, r
        return p1

    def poly_pow_mod(base, exp, mod_p):
        """Calculates (base^exp) mod mod_p for polynomials over F_m."""
        res = [1]
        base = normalize(list(base))
        mod_p = normalize(list(mod_p))
        
        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, base)
                _, res = poly_divmod(res, mod_p)
            base = poly_mul(base, base)
            _, base = poly_divmod(base, mod_p)
            exp //= 2
        return res

    A = []
    F = range(m)
    
    # Iterate through all possible values of 'a' in F
    for a in F:
        # P(x) = x^5 + ax + 3
        p = [3, a, 0, 0, 0, 1]

        # Step 1: Check for roots in F (test for factors of degree 1)
        has_root = False
        for c in F:
            val = 0
            # Evaluate P(c) using Horner's method
            for coeff in reversed(p):
                val = (val * c + coeff) % m
            if val == 0:
                has_root = True
                break
        if has_root:
            continue

        # Step 2: If no roots, check for irreducible factors of degree 2
        # An irreducible f(x) of degree > 2 has gcd(f(x), x^(q^2)-x) = 1
        x_poly = [0, 1]  # Represents the polynomial 'x'
        x_pow_q2 = poly_pow_mod(x_poly, m**2, p)
        rem = poly_sub(x_pow_q2, x_poly)
        g = poly_gcd(rem, p)
        
        # If gcd is a non-zero constant (degree 0, length 1), P(x) is irreducible
        if len(g) == 1 and g[0] != 0:
            A.append(a)
    
    # Step 3: Perform the final calculation
    max_A = max(A)
    min_A = min(A)
    len_A = len(A)
    result = max_A ** min_A - len_A
    
    print(f"The set of values 'a' for which the polynomial is irreducible is A = {sorted(A)}")
    print(f"The maximum element of A is max(A) = {max_A}")
    print(f"The minimum element of A is min(A) = {min_A}")
    print(f"The size of A is |A| = {len_A}")
    print("\nThe final calculation is max(A)^min(A) - |A|")
    print(f"The final equation is: {max_A}^{min_A} - {len_A} = {result}")

solve_polynomial_problem()