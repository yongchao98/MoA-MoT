def poly_mult(p1, p2, deg):
    """Multiplies two polynomials p1 and p2, returning the result up to a given degree."""
    n1 = len(p1)
    n2 = len(p2)
    res = [0] * (deg + 1)
    for i in range(n1):
        for j in range(n2):
            if i + j <= deg:
                res[i + j] += p1[i] * p2[j]
    return res

def poly_inv(p, deg):
    """Computes the multiplicative inverse of a polynomial p up to a given degree."""
    if p[0] == 0:
        raise ValueError("Cannot invert a polynomial with a zero constant term.")
    
    # Normalize p so that the constant term is 1
    const_term = p[0]
    p_norm = [c / const_term for c in p]
    
    res = [0] * (deg + 1)
    res[0] = 1.0 / const_term
    
    for i in range(1, deg + 1):
        s = 0
        for j in range(1, i + 1):
            if j < len(p_norm):
                s += p_norm[j] * res[i - j]
        res[i] = -s
        
    return [round(c) for c in res]

def solve():
    """
    Calculates the number of 1324-avoiding permutations of length n with k inversions for large n.
    This value, av_n^k(1324), is constant for n >= (4-1)*k.
    For k=3, this holds for n >= 9.
    The value is the coefficient of x^3 in the generating function:
    A(x) = (1/(1-x)) * product_{i>=1} (1 / (1 - x^i * C(x^i)))
    """
    DEGREE = 3

    # Catalan numbers: C(z) = 1 + z + 2z^2 + 5z^3 + ...
    C_poly = [1, 1, 2, 5]

    # Term 1: P1 = 1/(1-x)
    P1 = [1] * (DEGREE + 1)

    # Term 2: P2 = 1 / (1 - x*C(x))
    # x*C(x) = x + x^2 + 2x^3 + ...
    x_C_x_poly = [0] + C_poly[:-1]
    denom2 = [1] + [-c for c in x_C_x_poly[1:]]
    P2 = poly_inv(denom2, DEGREE)

    # Term 3: P3 = 1 / (1 - x^2*C(x^2))
    # C(x^2) = 1 + x^2 + ...
    # x^2*C(x^2) = x^2 + x^4 + ...
    x2_C_x2_poly = [0, 0, 1] # x^2 * (1)
    denom3 = [1, 0, -1]
    P3 = poly_inv(denom3, DEGREE)

    # Term 4: P4 = 1 / (1 - x^3*C(x^3))
    # C(x^3) = 1 + ...
    # x^3*C(x^3) = x^3 + ...
    x3_C_x3_poly = [0, 0, 0, 1] # x^3 * (1)
    denom4 = [1, 0, 0, -1]
    P4 = poly_inv(denom4, DEGREE)
    
    # A(x) = P1 * P2 * P3 * P4 * ... (mod x^4)
    temp_res = poly_mult(P1, P2, DEGREE)
    temp_res = poly_mult(temp_res, P3, DEGREE)
    final_poly = poly_mult(temp_res, P4, DEGREE)
    
    result = final_poly[DEGREE]
    
    # Print the equation representing the final answer
    print(f"av_{{333}}^{{3}}(1324) = {result}")

solve()