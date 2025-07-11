import numpy as np
import sympy
from sympy import Poly, Symbol

def get_bezout_matrix(p1, p2, gen):
    """
    Computes the Bezout matrix of two polynomials.
    This function calculates the coefficients of the Bezoutian (p1(x)p2(y) - p1(y)p2(x))/(x-y).
    """
    x = gen
    y = Symbol('y', real=True)

    # Convert Poly objects to expressions for substitution
    p1_expr = p1.as_expr()
    p2_expr = p2.as_expr()

    p1_y = p1_expr.subs(x, y)
    p2_y = p2_expr.subs(x, y)

    # Calculate the Bezoutian polynomial
    bezoutian = sympy.simplify((p1_expr * p2_y - p1_y * p2_expr) / (x - y))
    
    n = max(p1.degree(), p2.degree())
    B = sympy.zeros(n, n)
    
    # Extract coefficients to form the Bezout matrix
    poly_bezoutian = Poly(bezoutian, x, y)
    for i in range(n):
        for j in range(n):
            B[i, j] = poly_bezoutian.coeff_monomial(x**i * y**j)
            
    return B

def get_sylvester_matrix(p1, p2, gen):
    """
    Computes the Sylvester matrix using SymPy's built-in function.
    """
    return sympy.sylvester(p1, p2, gen=gen)

def main():
    """
    Main function to execute the full calculation pipeline.
    """
    # Step 1: Determine the coefficients C_i based on interpretation of visualizations
    
    # For V_3 (interpreted as representing a group of order 4, C2xC2):
    # R_j3 = [Sum Char Table, # irreps, Order, Exponent] = [4, 4, 4, 2]
    # Sum(R_j3) = 14, Sum(R_j3^2) = 16+16+16+4 = 52
    C3 = int(np.floor(52 / 14))

    # For 10-dimensional objects (V_1, V_5, V_6, V_7, V_9), assuming group C_10:
    # R_ji = [10, 10, 10, 10]
    # Sum(R_ji) = 40, Sum(R_ji^2) = 100+100+100+100 = 400
    Ci_10 = int(np.floor(400 / 40))
    C1 = C5 = C6 = C7 = C9 = Ci_10

    # For 20-dimensional objects (V_2, V_4, V_8), assuming group C_20:
    # R_ji = [20, 20, 20, 20]
    # Sum(R_ji) = 80, Sum(R_ji^2) = 400+400+400+400 = 1600
    Ci_20 = int(np.floor(1600 / 80))
    C2 = C4 = C8 = Ci_20
    
    C = [C1, C2, C3, C4, C5, C6, C7, C8, C9]

    # Step 2: Construct polynomials P(x), Q(x), and S(x)
    x = Symbol('x')
    P_x = sum(c * x**i for i, c in enumerate(C, 1))

    i = sympy.I
    P_ix = P_x.subs(x, i * x).expand()
    
    Q_x = sympy.re(P_ix)
    S_x = sympy.im(P_ix)

    Q_poly = Poly(Q_x, x)
    S_poly = Poly(S_x, x)
    
    # Step 3: Construct Bezout and Sylvester matrices
    
    # M1 = Bm[Q(x), S(x), x]
    M1 = get_bezout_matrix(Q_poly, S_poly, x)

    # M2 = Sm[Q(x), x^10 + S(x), x]
    S2_poly = Poly(x**10 + S_x, x)
    M2 = get_sylvester_matrix(Q_poly, S2_poly, x)

    # Step 4: Compute the final trace T
    trace_M1 = M1.trace()
    trace_M2 = M2.trace()
    
    T = trace_M1 * 2 + trace_M2
    
    # Step 5: Print the results in a structured format
    q_poly_coeffs = Q_poly.all_coeffs()
    s_poly_coeffs = S_poly.all_coeffs()
    
    q_str_parts = []
    for j, c in enumerate(q_poly_coeffs):
        if c != 0: q_str_parts.append(f"{int(c)}*x^{Q_poly.degree() - j}")
    q_str = " + ".join(q_str_parts).replace('+ -', '- ')

    s_str_parts = []
    for j, c in enumerate(s_poly_coeffs):
        if c != 0: s_str_parts.append(f"{int(c)}*x^{S_poly.degree() - j}")
    s_str = " + ".join(s_str_parts).replace('+ -', '- ')


    print("Step 1: Coefficients C_i")
    print(f"C = {C}")
    print("\nStep 2: Polynomials")
    print(f"Q(x) = {q_str}")
    print(f"S(x) = {s_str}")
    print("\nStep 3: Bezout Matrix M1 Trace")
    print(f"Tr(M1) = {trace_M1}")
    print("\nStep 4: Sylvester Matrix M2 Trace")
    print(f"Tr(M2) = {trace_M2}")
    print("\nStep 5: Final Equation for T")
    print("T = Tr(M1 \u2297 I2 + M2) = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * ({trace_M1}) + ({trace_M2})")
    print(f"T = {2*trace_M1} + ({trace_M2})")
    print(f"T = {T}")
    print("\nFinal Answer:")
    print(f"<<<{T}>>>")

if __name__ == "__main__":
    main()
