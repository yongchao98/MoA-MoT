import sympy as sp

def solve():
    """
    Solves the problem by calculating coefficients, defining polynomials,
    and computing the trace of the final matrix expression.
    """
    # Step 1: Calculate coefficients C_i based on problem analysis.
    # C_i is the floor of the contraharmonic mean of R_ji values.
    # From visual analysis, V3 has order 4, others have order 10.
    # This difference is the source of the non-uniformity in C_i values.
    
    C = {}
    # For i != 3, any combination of determinable properties (value 10)
    # and non-determinable ones (value 0) leads to a contraharmonic mean of 10.
    for i in [1, 2, 4, 5, 6, 7, 8, 9]:
        C[i] = 10
    
    # For i = 3, we extract R3 = 4, and other R values are 0.
    # Contraharmonic mean = (4^2) / 4 = 4.
    C[3] = 4

    print("Calculated coefficients C_i:")
    for i in range(1, 10):
        print(f"C_{i} = {C[i]}")

    # Step 2: Construct polynomials P(x), Q(x), S(x)
    x = sp.Symbol('x')
    P_x = sum(C[i] * x**i for i in range(1, 10))

    # P(ix) = Q(x) + i*S(x)
    P_ix = P_x.subs(x, sp.I * x)
    Q_x = sp.re(P_ix).expand()
    S_x = sp.im(P_ix).expand()

    print(f"\nPolynomial P(x) = {P_x}")
    print(f"Real part of P(ix), Q(x) = {Q_x}")
    print(f"Imaginary part of P(ix), S(x) = {S_x}")

    # Step 3: Calculate traces using mathematical theorems.
    # M1 = Bm[Q(x), S(x), x], M2 = Sm[Q(x), x^10 + S(x), x]
    
    # Theorem for trace of Bezout matrix Bm(p,q):
    # If deg(p) < deg(q) = n, then Tr(Bm(p,q)) is the coefficient of x^(n-1) in p(x).
    p_bez = Q_x # deg 8
    q_bez = S_x # deg 9, so n = 9
    tr_M1 = p_bez.coeff(x, q_bez.degree(x) - 1)

    # Theorem for trace of Sylvester matrix Sm(p,q):
    # Tr(Sm(p,q)) = n*a_m + m*b_0 where p has deg m, lead coeff a_m
    # and q has deg n, const term b_0.
    p_syl = Q_x
    q_syl = x**10 + S_x
    m_syl = p_syl.degree(x)
    n_syl = q_syl.degree(x)
    a_m = p_syl.coeff(x, m_syl)
    b_0 = q_syl.coeff(x, 0)
    tr_M2 = n_syl * a_m + m_syl * b_0

    # Step 4: Final trace calculation
    # T = Tr(M1 x I2 + M2) = 2 * Tr(M1) + Tr(M2)
    
    print("\n--- Final Calculation ---")
    print(f"Trace of Bezout Matrix M1: Tr(M1) = {tr_M1}")
    print(f"Trace of Sylvester Matrix M2: Tr(M2) = {tr_M2}")
    print(f"The final equation is: T = 2 * Tr(M1) + Tr(M2)")
    print(f"Substituting the values: T = 2 * {tr_M1} + {tr_M2}")
    final_T = 2 * tr_M1 + tr_M2
    print(f"The final trace is: T = {final_T}")
    print(f"\n<<<120>>>")

solve()