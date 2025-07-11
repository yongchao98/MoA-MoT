import math

def solve_polynomial_problem():
    """
    Solves the problem by finding irreducible polynomials and performing the final calculation.
    """
    # The finite field is F_7, so all arithmetic is modulo 7.
    F = 7

    # --- Helper functions for polynomial arithmetic over F ---

    def poly_norm(p):
        """Removes leading zeros from a polynomial represented as a list."""
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p

    def poly_add(p1, p2):
        """Adds two polynomials p1, p2 over F."""
        n = max(len(p1), len(p2))
        res = [0] * n
        for i in range(len(p1)):
            res[i] = p1[i]
        for i in range(len(p2)):
            res[i] = (res[i] + p2[i]) % F
        return poly_norm(res)

    def poly_mul(p1, p2):
        """Multiplies two polynomials p1, p2 over F."""
        n1, n2 = len(p1), len(p2)
        res = [0] * (n1 + n2 - 1)
        for i in range(n1):
            for j in range(n2):
                res[i+j] = (res[i+j] + p1[i] * p2[j]) % F
        return poly_norm(res)

    def mod_inverse(n):
        """Computes the modular inverse of n modulo F."""
        return pow(n, F - 2, F)

    def poly_div(N, D):
        """
        Performs polynomial long division N / D over F.
        Returns (quotient, remainder).
        """
        N = list(N)
        D = list(D)
        if D == [0]:
            raise ZeroDivisionError("Polynomial division by zero")
        
        degN, degD = len(N) - 1, len(D) - 1
        if degN < degD:
            return [0], N

        q = [0] * (degN - degD + 1)
        inv_d_lead = mod_inverse(D[-1])

        while degN >= degD:
            d = [0] * (degN - degD) + D
            mult = (N[-1] * inv_d_lead) % F
            q[degN - degD] = mult
            
            d_mult = [(mult * c) % F for c in d]
            
            N = poly_add(N, [(F - c) % F for c in d_mult])
            degN = len(N) - 1
            
        return poly_norm(q), N


    def poly_gcd(p1, p2):
        """Computes the GCD of two polynomials over F using the Euclidean algorithm."""
        while p2 != [0]:
            _, r = poly_div(p1, p2)
            p1, p2 = p2, r
        
        # Normalize to make the GCD monic
        if p1 and p1[-1] != 0:
            inv_lead = mod_inverse(p1[-1])
            p1 = [(c * inv_lead) % F for c in p1]
        return p1

    def poly_pow_mod(base, exp, mod_poly):
        """Computes (base^exp) % mod_poly over F using binary exponentiation."""
        res = [1]
        base = list(base)
        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, base)
                _, res = poly_div(res, mod_poly)
            base = poly_mul(base, base)
            _, base = poly_div(base, mod_poly)
            exp //= 2
        return res

    def is_irreducible(p):
        """
        Checks if polynomial p is irreducible over F using the GCD test.
        A polynomial f(x) of degree n is irreducible iff
        gcd(f(x), x^(q^d) - x) = 1 for all 1 <= d <= n/2.
        """
        deg_p = len(p) - 1
        if deg_p <= 1:
            return True if deg_p > 0 else False
            
        x_poly = [0, 1]  # The polynomial x

        for d in range(1, deg_p // 2 + 1):
            q_power_d = F**d
            
            # Compute x^(F^d) mod p
            term = poly_pow_mod(x_poly, q_power_d, p)
            
            # Then compute h = x^(F^d) - x (mod p)
            h = poly_add(term, [0, F - 1])
            
            g = poly_gcd(p, h)

            # If gcd is not 1 (a constant polynomial), it's reducible
            if g != [1]:
                return False
                
        return True

    # --- Main Logic ---
    A = []
    
    # Iterate through all possible values of 'a' in the field F_7
    for a in range(F):
        # The polynomial is p(x) = x^5 + ax + 3
        # In list form (coefficients of x^0, x^1, ...): [3, a, 0, 0, 0, 1]
        p = [3, a, 0, 0, 0, 1]
        
        if is_irreducible(p):
            A.append(a)

    print(f"The finite field is F = Z_7.")
    print(f"The polynomial is p(x) = x^5 + ax + 3.")
    print(f"The set A of 'a' values for which p(x) is irreducible is: {A}")

    if not A:
        print("The set A is empty, so the calculation cannot be performed.")
        return

    min_A = min(A)
    max_A = max(A)
    cardinality_A = len(A)

    result = max_A ** min_A - cardinality_A

    print(f"\nThe minimum element of A is min(A) = {min_A}")
    print(f"The maximum element of A is max(A) = {max_A}")
    print(f"The cardinality of A is |A| = {cardinality_A}")
    print("\nThe required calculation is max(A)^min(A) - |A|.")
    print(f"The result is {max_A}^{min_A} - {cardinality_A} = {result}")

solve_polynomial_problem()
<<<3>>>