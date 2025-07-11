import sys

# Set recursion limit higher for polynomial operations, though likely not necessary
# for this specific problem size.
sys.setrecursionlimit(2000)

def solve():
    """
    Solves the problem by finding the set A and computing the final expression.
    """
    F_SIZE = 7

    # Utility functions for polynomial arithmetic over F_SIZE.
    # Polynomials are represented as lists of coefficients, e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2.

    def poly_trim(p):
        """Removes leading zero coefficients."""
        while len(p) > 1 and p[-1] == 0:
            p.pop()
        return p

    def poly_add(p, q):
        """Adds two polynomials."""
        n = max(len(p), len(q))
        res = [0] * n
        for i in range(n):
            c1 = p[i] if i < len(p) else 0
            c2 = q[i] if i < len(q) else 0
            res[i] = (c1 + c2) % F_SIZE
        return poly_trim(res)

    def poly_sub(p, q):
        """Subtracts polynomial q from p."""
        n = max(len(p), len(q))
        res = [0] * n
        for i in range(n):
            c1 = p[i] if i < len(p) else 0
            c2 = q[i] if i < len(q) else 0
            res[i] = (c1 - c2 + F_SIZE) % F_SIZE
        return poly_trim(res)

    def poly_mul(p, q):
        """Multiplies two polynomials."""
        if (len(p) == 1 and p[0] == 0) or (len(q) == 1 and q[0] == 0):
            return [0]
        n, m = len(p), len(q)
        res = [0] * (n + m - 1)
        for i in range(n):
            for j in range(m):
                res[i + j] = (res[i + j] + p[i] * q[j]) % F_SIZE
        return poly_trim(res)

    def poly_long_div(dividend, divisor):
        """Performs long division, returns (quotient, remainder)."""
        dividend, divisor = list(dividend), list(divisor)
        poly_trim(dividend)
        poly_trim(divisor)

        if len(divisor) == 1 and divisor[0] == 0:
            raise ZeroDivisionError("Polynomial division by zero")

        if len(dividend) < len(divisor):
            return ([0], dividend)

        inv_lead_d = pow(divisor[-1], F_SIZE - 2, F_SIZE)
        deg_n, deg_d = len(dividend) - 1, len(divisor) - 1
        quotient = [0] * (deg_n - deg_d + 1)
        
        while len(dividend) >= len(divisor):
            deg_rem = len(dividend) - 1
            lead_rem_coeff = dividend[-1]
            
            q_deg = deg_rem - deg_d
            q_coeff = (lead_rem_coeff * inv_lead_d) % F_SIZE
            quotient[q_deg] = q_coeff
            
            term_to_sub = [0] * (q_deg + 1)
            term_to_sub[q_deg] = q_coeff
            
            subtrahend = poly_mul(divisor, term_to_sub)
            dividend = poly_sub(dividend, subtrahend)
            poly_trim(dividend)
            
        return (quotient, dividend)

    def poly_gcd(p, q):
        """Computes GCD of two polynomials using the Euclidean algorithm."""
        while not (len(q) == 1 and q[0] == 0):
            p, q = q, poly_long_div(p, q)[1]
        return p

    def poly_pow_mod(base, exp, modulus):
        """Computes (base^exp) mod modulus for polynomials."""
        res = [1]
        while exp > 0:
            if exp % 2 == 1:
                res = poly_mul(res, base)
                res = poly_long_div(res, modulus)[1]
            base = poly_mul(base, base)
            base = poly_long_div(base, modulus)[1]
            exp //= 2
        return res

    def evaluate_poly(p, x):
        """Evaluates polynomial p at value x using Horner's method."""
        res = 0
        for i in range(len(p) - 1, -1, -1):
            res = (res * x + p[i]) % F_SIZE
        return res

    A = []
    # Iterate through all possible values of 'a' in F_7
    for a in range(F_SIZE):
        # p(x) = x^5 + ax + 3
        p = [3, a, 0, 0, 0, 1]

        # Test 1: Check for roots in F_7. If it has roots, it's reducible.
        has_root = False
        for x_val in range(F_SIZE):
            if evaluate_poly(p, x_val) == 0:
                has_root = True
                break
        if has_root:
            continue

        # Test 2: Check for irreducible quadratic factors.
        # This is true if gcd(p(x), x^(7^2) - x) is a constant.
        x_poly = [0, 1]
        h_poly = poly_pow_mod(x_poly, F_SIZE**2, p)
        g_poly = poly_sub(h_poly, x_poly)
        gcd_res = poly_gcd(p, g_poly)
        
        # A polynomial is irreducible if the gcd is a non-zero constant (degree 0).
        if len(gcd_res) == 1 and gcd_res[0] != 0:
            A.append(a)

    # Final computation
    if not A:
        print("The set A is empty.")
        return
        
    min_A = min(A)
    max_A = max(A)
    size_A = len(A)
    
    result = max_A**min_A - size_A
    
    print(f"Let F be the finite field of order 7.")
    print(f"Let A be the set of elements 'a' in F such that the polynomial x^5+ax+3 is irreducible.")
    print("")
    print(f"The set A is: {A}")
    print(f"The minimum of A is: min(A) = {min_A}")
    print(f"The maximum of A is: max(A) = {max_A}")
    print(f"The size of A is: |A| = {size_A}")
    print("")
    print(f"The required calculation is max(A)^min(A) - |A|:")
    print(f"= {max_A}^{min_A} - {size_A}")
    print(f"= {max_A**min_A} - {size_A}")
    print(f"= {result}")
    
solve()
<<<3>>>