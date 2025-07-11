def poly_deg(p):
    """
    Calculates the degree of a polynomial represented as a list of coefficients.
    The zero polynomial has degree -1.
    """
    i = len(p) - 1
    while i >= 0 and p[i] == 0:
        i -= 1
    return i

def poly_divmod(num, den, m):
    """
    Performs polynomial long division over a finite field F_m.
    Args:
        num (list): The numerator polynomial.
        den (list): The denominator polynomial.
        m (int): The modulus of the field.
    Returns:
        tuple: A tuple containing the quotient and remainder polynomials.
    """
    # Make copies to avoid modifying original lists
    num = list(num)
    den = list(den)

    deg_den = poly_deg(den)
    if deg_den == -1:
        raise ZeroDivisionError("Polynomial division by zero")

    deg_num = poly_deg(num)
    
    if deg_num < deg_den:
        return ([0], num)

    # Calculate the modular inverse of the leading coefficient of the denominator
    inv_lc_den = pow(den[deg_den], m - 2, m)

    q = [0] * (deg_num - deg_den + 1)
    
    for i in range(deg_num, deg_den - 1, -1):
        if num[i] != 0:
            coeff = (num[i] * inv_lc_den) % m
            q_deg = i - deg_den
            q[q_deg] = coeff

            for j in range(deg_den + 1):
                num[i - deg_den + j] = (num[i - deg_den + j] - coeff * den[j]) % m
    
    r = num
    return q, r

def poly_gcd(a, b, m):
    """
    Calculates the greatest common divisor of two polynomials over F_m
    using the Euclidean algorithm.
    """
    while poly_deg(b) != -1:
        a, b = b, poly_divmod(a, b, m)[1]
    return a

def solve():
    """
    Solves the problem by finding the set A and performing the final calculation.
    """
    F_ord = 7
    A = []

    # The polynomial is x^(7^2) - x = x^49 - x
    poly_x49_minus_x = [0] * 50
    poly_x49_minus_x[49] = 1
    poly_x49_minus_x[1] = F_ord - 1

    # Iterate through all possible values for 'a' in F_7
    for a in range(F_ord):
        # Construct the polynomial P_a(x) = x^5 + ax + 3
        p_a = [0] * 6
        p_a[0] = 3
        p_a[1] = a
        p_a[5] = 1

        # A polynomial of degree 5 is irreducible over F_7 iff
        # gcd(P_a(x), x^(7^2) - x) is a constant (degree 0).
        gcd_res = poly_gcd(p_a, poly_x49_minus_x, F_ord)
        
        if poly_deg(gcd_res) == 0:
            A.append(a)

    # Perform the final calculation: max(A)^min(A) - |A|
    min_A = min(A)
    max_A = max(A)
    len_A = len(A)
    
    result = max_A ** min_A - len_A
    
    # Print the required output
    print(f"The set A of coefficients 'a' for which the polynomial is irreducible is: {A}")
    print(f"The calculation is max(A)^min(A) - |A|")
    print(f"Substituting the values: {max_A}^{min_A} - {len_A} = {result}")

solve()
<<<3>>>