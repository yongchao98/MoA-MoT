import numpy as np
import sympy as sp

def calculate_final_trace():
    """
    This function performs the mathematical computations as outlined in the problem.
    """
    # Step 2: Define the properties of the identified group D_5
    R_values = [8, 4, 10, 10]  # Sum of chars, #classes, order, exponent

    # Step 3: Compute C_i
    r_sq = [r**2 for r in R_values]
    contraharmonic_mean = sum(r_sq) / sum(R_values)
    C = int(np.floor(contraharmonic_mean))

    # Step 4: Construct polynomials Q(x) and S(x)
    x = sp.Symbol('x')
    P = sum(C * x**i for i in range(1, 10))
    Pix = P.subs(x, sp.I * x)
    
    # Sympy can separate the real (Q) and imaginary (S) parts
    Qx_expr = sp.re(Pix).expand()
    Sx_expr = sp.im(Pix.as_expr()).expand() # im(poly) is tricky, so use im(expr)
    
    Qx = sp.Poly(Qx_expr, x)
    Sx = sp.Poly(Sx_expr, x)

    # Step 5: Calculate traces of M1 and M2
    
    # Trace of M2 = Sm(Q, x^10 + S)
    # A = Q(x), n = deg(Q) = 8, leading coeff a_n = 8
    # B = x^10 + S(x), m = deg(B) = 10, constant term b_0 = 0
    n = Qx.degree()
    a_n = Qx.coeffs()[0]
    
    Bx = sp.Poly(x**10 + Sx_expr, x)
    m = Bx.degree()
    b_0 = Bx.coeffs()[-1] if Bx.coeffs()[-1].is_number else 0
    
    tr_m2 = m * a_n + n * b_0
    
    # Trace of M1 = Bm(Q, S)
    # We use a formula for trace of Bezoutian for odd/even polynomials.
    # Tr(M1) = q2(s3-s1) + q4(s5-s3) + q6(s7-s5) + q8(s9-s7)
    q_coeffs = Qx.all_coeffs() # in descending power order
    s_coeffs = Sx.all_coeffs()
    q = lambda i: q_coeffs[Qx.degree()-i] if i <= Qx.degree() and i >= 0 else 0
    s = lambda i: s_coeffs[Sx.degree()-i] if i <= Sx.degree() and i >= 0 else 0
    
    tr_m1 = q(2)*(s(3)-s(1)) + q(4)*(s(5)-s(3)) + q(6)*(s(7)-s(5)) + q(8)*(s(9)-s(7))

    # Step 6: Final Calculation
    T = 2 * tr_m1 + tr_m2
    
    print(f"R values: {R_values}")
    print(f"C: floor({contraharmonic_mean:.2f}) = {C}")
    print(f"Q(x) = {Qx_expr}")
    print(f"S(x) = {Sx_expr}")
    print(f"Tr(M1) = {tr_m1}")
    print(f"Tr(M2) = {tr_m2}")
    print(f"Final equation: 2 * {tr_m1} + {tr_m2} = {T}")
    print(f"The trace T is {T}")

calculate_final_trace()