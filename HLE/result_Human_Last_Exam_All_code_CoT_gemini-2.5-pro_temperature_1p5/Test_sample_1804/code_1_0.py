import math

def poly_eval(poly, x, p):
    """Evaluates a polynomial at x in a field of order p."""
    res = 0
    for coeff in poly:
        res = (res * x + coeff) % p
    return res

def poly_divmod(num, den, p):
    """
    Performs polynomial long division in a field of order p.
    Returns quotient and remainder.
    Coefficients are given from highest degree to lowest.
    """
    num = list(num)
    den = list(den)
    
    if not den:
        raise ZeroDivisionError("Divisor cannot be empty")
    
    deg_num = len(num) - 1
    deg_den = len(den) - 1

    if deg_num < deg_den:
        return [0], num

    quotient = [0] * (deg_num - deg_den + 1)
    
    # Modular inverse of leading coefficient of the divisor
    try:
        inv_lead_den = pow(den[0], -1, p)
    except ValueError:
        raise ZeroDivisionError("Leading coefficient of divisor has no modular inverse.")

    while deg_num >= deg_den:
        d = [0] * (deg_num - deg_den) + den
        
        # Calculate multiplier for current subtraction step
        mult = (num[0] * inv_lead_den) % p
        quotient[len(num) - len(d)] = mult

        # Subtract d*mult from num
        for i in range(len(d)):
            num[i] = (num[i] - d[i] * mult) % p
        
        # Remove leading zero of the remainder
        num.pop(0)
        deg_num -= 1
        
    # Clean up leading zeros in remainder
    while len(num) > 1 and num[0] == 0:
        num.pop(0)
    
    # If remainder is empty, it's zero
    if not num:
        num = [0]
        
    return quotient, num

def solve_finite_field_problem():
    """
    Solves the problem described by the user.
    """
    p = 7
    F = list(range(p))

    # Step 1: Find all monic irreducible quadratic polynomials over F_7
    # q(x) = x^2 + bx + c is irreducible if it has no roots in F_7.
    irred_quadratics = []
    for b in F:
        for c in F:
            q = [1, b, c]
            has_root = False
            for x in F:
                if poly_eval(q, x, p) == 0:
                    has_root = True
                    break
            if not has_root:
                irred_quadratics.append(q)

    # Step 2: Find the set A for the polynomial x^5 + ax + 3
    set_A = []
    for a in F:
        poly_p = [1, 0, 0, 0, a, 3]
        
        # 2a: Check for roots in F_7
        is_reducible_by_root = False
        for x in F:
            if poly_eval(poly_p, x, p) == 0:
                is_reducible_by_root = True
                break
        
        if is_reducible_by_root:
            continue
            
        # 2b: If no roots, check for irreducible quadratic factors
        is_reducible_by_factor = False
        for q in irred_quadratics:
            _quotient, remainder = poly_divmod(poly_p, q, p)
            if remainder == [0]:
                is_reducible_by_factor = True
                break
        
        if is_reducible_by_factor:
            continue

        # If no roots and no irreducible quadratic factors, it's irreducible
        set_A.append(a)

    # Step 3: Calculate the final expression
    if not set_A:
        print("Set A is empty, cannot perform calculation.")
        return
        
    max_A = max(set_A)
    min_A = min(set_A)
    card_A = len(set_A)

    result = int(math.pow(max_A, min_A)) - card_A

    print(f"The set A is {set_A}.")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A| = {card_A}")
    print(f"The result of max(A)^min(A) - |A| is:")
    print(f"{max_A}^{min_A} - {card_A} = {result}")

solve_finite_field_problem()