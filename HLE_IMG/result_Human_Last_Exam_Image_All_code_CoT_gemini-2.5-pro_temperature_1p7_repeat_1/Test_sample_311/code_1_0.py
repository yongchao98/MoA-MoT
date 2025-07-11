import numpy as np
from sympy import symbols, Poly, expand, Matrix

def solve_problem():
    """
    Solves the problem by calculating the final trace T based on derived coefficients.
    """
    # Step 1: Define the coefficients C_i derived from the analysis.
    # C_i for i=1 to 9
    C = [0, 8, 5, 3, 16, 8, 10, 8, 10, 10]

    # Step 2: Construct polynomials P(x), Q(x), S(x)
    x, y = symbols('x y')
    P_coeffs_rev = list(reversed(C))
    P = Poly(P_coeffs_rev, x)

    # Substitute x with i*x and find real (Q) and imaginary (S) parts
    Pix_expanded = expand(P.subs(x, 1j * x))
    Q_poly = Poly(Pix_expanded.as_real_imag()[0], x)
    S_poly = Poly(Pix_expanded.as_real_imag()[1], x)

    # Step 3: Calculate Trace of M1 (Bezout Matrix)
    # The Bezoutian B(Q,S) is a 9x9 matrix defined by:
    # (Q(x)S(y) - Q(y)S(x)) / (x - y) = sum_{i,j=0 to 8} b_ij * x^i * y^j
    # We can compute the trace manually from polynomial coefficients.
    q_coeffs = Q_poly.all_coeffs() # High-to-low degree
    s_coeffs = S_poly.all_coeffs() # High-to-low degree
    
    # Pad with zeros to get standard coefficient lists (power 0 to n-1)
    n = max(Q_poly.degree(), S_poly.degree())
    q = np.zeros(n + 1)
    s = np.zeros(n + 1)
    
    # sympy coeffs are reversed (high power first)
    for i, coeff in enumerate(reversed(Q_poly.all_coeffs())):
        q[i] = coeff
    for i, coeff in enumerate(reversed(S_poly.all_coeffs())):
        s[i] = coeff

    # Calculate diagonal entries of Bezout matrix
    trace_m1 = 0
    for i in range(n):
        b_ii = 0
        for k in range(i + 1):
            s_idx = i - k
            q_idx1 = i + k + 1
            q_idx2 = i - k
            s_idx2 = i + k + 1
            term1 = s[s_idx] * q[q_idx1] if s_idx < len(s) and q_idx1 < len(q) else 0
            term2 = q[q_idx2] * s[s_idx2] if q_idx2 < len(q) and s_idx2 < len(s) else 0
            b_ii += term1 - term2
        trace_m1 += b_ii
    trace_m1 = int(trace_m1)
    
    # Step 4: Calculate Trace of M2 (Sylvester Matrix)
    # Tr(Syl(p,q)) is the sum of n copies of the leading coefficient of p
    # and m copies of the constant term of q.
    p = Q_poly
    q_sylvester = Poly(x**10) + S_poly
    
    m_p = p.degree()
    n_q = q_sylvester.degree()
    
    lead_coeff_p = p.LC()
    const_term_q = q_sylvester.eval(0)
    
    trace_m2 = n_q * lead_coeff_p + m_p * const_term_q
    trace_m2 = int(trace_m2)

    # Step 5: Final calculation of T
    T = 2 * trace_m1 + trace_m2

    print(f"Coefficients C_i: {[C[i] for i in range(1, 10)]}")
    print(f"Q(x) = {Q_poly.as_expr()}")
    print(f"S(x) = {S_poly.as_expr()}")
    print(f"Trace(M1) = {trace_m1}")
    print(f"Trace(M2) = {trace_m2}")
    print(f"T = Tr(M1 x I2 + M2)")
    print(f"T = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * ({trace_m1}) + ({trace_m2})")
    print(f"T = {2 * trace_m1} + {trace_m2}")
    print(f"T = {T}")
    print(f"<<<{T}>>>")

solve_problem()