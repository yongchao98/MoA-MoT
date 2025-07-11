import sympy
import numpy as np

def solve_problem():
    """
    Solves the multi-step problem based on the properties of the D4 group.
    """
    # Step 1: Define group properties for D4.
    R1 = 8  # Sum of character table entries
    R2 = 5  # Number of irreducible representations
    R3 = 8  # Order of the group
    R4 = 4  # Exponent of the group

    # Step 2: Calculate the contraharmonic mean constant C.
    R_vals = np.array([R1, R2, R3, R4])
    sum_of_squares = np.sum(R_vals**2)
    sum_of_vals = np.sum(R_vals)
    contraharmonic_mean = sum_of_squares / sum_of_vals
    C = int(np.floor(contraharmonic_mean))

    # Step 3: Construct polynomials P(x), Q(x), and S(x).
    x = sympy.Symbol('x')
    P_x = sum(C * x**i for i in range(1, 10))
    P_ix = P_x.subs(x, sympy.I * x)
    Q_x = sympy.re(P_ix)
    S_x = sympy.im(P_ix)

    # Step 4: Compute matrix traces.
    
    # Trace of M1 = Bm[Q(x), S(x), x]
    Q_poly = sympy.Poly(Q_x, x)
    S_poly = sympy.Poly(S_x, x)
    
    # Helper function to compute the Bezout matrix for two polynomials.
    def bezout_matrix(p_poly, q_poly, var):
        k = max(p_poly.degree(), q_poly.degree())
        y = sympy.Symbol('y')
        num = p_poly.subs({var: x}) * q_poly.subs({var: y}) - p_poly.subs({var: y}) * q_poly.subs({var: x})
        # Use cancel to simplify the rational expression before division
        bezoutian_poly = sympy.poly(sympy.cancel(num / (x - y)), x, y)
        mat = sympy.zeros(k, k)
        for i in range(k):
            for j in range(k):
                mat[i, j] = bezoutian_poly.coeff_monomial(x**i * y**j)
        return mat

    M1 = bezout_matrix(Q_poly, S_poly, x)
    tr_M1 = M1.trace()

    # Trace of M2 = Sm[Q(x), x^10 + S(x), x]
    P1_poly = Q_poly
    P2_poly = sympy.Poly(x**10 + S_x, x)
    n = P1_poly.degree()
    m = P2_poly.degree()
    # Using the formula Tr(Sm(p,q)) = m*LC(p) + n*const(q)
    tr_M2 = m * P1_poly.LC() + n * P2_poly.coeff_monomial(x**0)

    # Step 5: Calculate the final trace T.
    T = 2 * tr_M1 + tr_M2
    
    # Output the final equation with all numbers
    print(f"The final trace T is computed from the equation: T = 2 * Tr(M1) + Tr(M2)")
    print(f"Tr(M1) = {tr_M1}")
    print(f"Tr(M2) = {tr_M2}")
    print(f"T = 2 * {tr_M1} + {tr_M2} = {T}")
    
    print(f"\nFinal Answer:")
    print(f"<<<{T}>>>")

solve_problem()