import sympy

def solve_crystallographic_problem():
    """
    This function implements the described method to solve the problem.
    """
    # Step 1 & 2: Determine group properties and calculate C.
    # Based on the analysis, the group is D5.
    # R1: Sum of char table entries = 8
    # R2: Number of irreps = 4
    # R3: Group order = 10
    # R4: Group exponent = 10
    R_values = [8, 4, 10, 10]
    
    sum_R = sum(R_values)
    sum_R_sq = sum(val**2 for val in R_values)
    contraharmonic_mean = sum_R_sq / sum_R
    C = int(sympy.floor(contraharmonic_mean))
    
    # Step 3: Construct polynomials
    x = sympy.Symbol('x')
    P_x = sum(C * x**i for i in range(1, 10))
    P_ix = P_x.subs(x, sympy.I * x)
    
    # Q(x) is the real part of P(ix), S(x) is the imaginary part.
    Q_x = sympy.re(P_ix)
    S_x = sympy.im(P_ix)

    # Step 4: Construct matrices M1 and M2
    Q_poly = sympy.Poly(Q_x, x)
    S_poly = sympy.Poly(S_x, x)
    
    # M1 is the Bezout matrix of S(x) and Q(x).
    # Per prompt, Bm[Q, S], but Sympy requires deg(f) >= deg(g) for bezout(f, g)
    # The trace is identical for Bm(Q,S) and Bm(S,Q).
    if S_poly.degree() < Q_poly.degree():
        S_poly, Q_poly = Q_poly, S_poly # Ensure first poly has higher degree
        
    M1 = sympy.polys.matrices.bezout(S_poly, Q_poly)
    
    # M2 is the Sylvester matrix of Q(x) and x^10 + S(x).
    T_x = x**10 + S_x
    T_poly = sympy.Poly(T_x, x)
    M2 = sympy.polys.matrices.sylvester(sympy.Poly(Q_x, x), T_poly)

    # Step 5: Calculate traces
    tr_M1 = M1.trace()
    tr_M2 = M2.trace()

    # Step 6: Final calculation of T
    # T = Tr(M1 (kronecker) I_2 + M2) = 2*Tr(M1) + Tr(M2)
    final_trace = 2 * tr_M1 + tr_M2

    # Print the breakdown of the final equation as requested
    print("--- Calculation Breakdown ---")
    print(f"R values (R1, R2, R3, R4): {R_values}")
    print(f"Contraharmonic Mean: {contraharmonic_mean:.4f}")
    print(f"Coefficient C = floor({contraharmonic_mean:.4f}) = {C}")
    print(f"Tr(M1) = {tr_M1}")
    print(f"Tr(M2) = {tr_M2}")
    print("\n--- Final Equation ---")
    print(f"T = 2 * Tr(M1) + Tr(M2)")
    print(f"T = 2 * ({tr_M1}) + ({tr_M2})")
    print(f"T = {2 * tr_M1} + {tr_M2}")
    print(f"T = {final_trace}")

solve_crystallographic_problem()
