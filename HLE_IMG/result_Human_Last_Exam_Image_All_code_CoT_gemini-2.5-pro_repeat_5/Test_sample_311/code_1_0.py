import sympy as sp

def solve():
    """
    This function implements the entire logic to solve the problem.
    It calculates the properties of the D5 group, computes the polynomial coefficients,
    constructs the necessary polynomials, and calculates the traces of the
    Bezout and Sylvester matrices to find the final answer T.
    """

    # Step 1: Define the properties of the identified group, D5.
    # R1: Sum of all entries within the Character Table.
    # For D5, this is 8.
    R1 = 8
    # R2: Number of irreducible representations.
    # For D5, this is 4.
    R2 = 4
    # R3: Order of the corresponding Crystallographic Point Group.
    # For D5, this is 10.
    R3 = 10
    # R4: Exponent of the corresponding Crystallographic Point Group.
    # For D5, this is lcm(orders of elements) = lcm(1,2,5) = 10.
    R4 = 10

    # Step 2: Compute the constant coefficient C.
    # C is the floor of the contraharmonic mean of the four quantities.
    R_values = [R1, R2, R3, R4]
    sum_R = sum(R_values)
    sum_sq_R = sum(r**2 for r in R_values)
    contraharmonic_mean = sum_sq_R / sum_R
    C = int(sp.floor(contraharmonic_mean))

    # Step 3: Construct the polynomial P(x) and find Q(x) and S(x).
    x, y = sp.symbols('x y')
    P_poly = sum(C * x**i for i in range(1, 10))
    
    # Find real (Q) and imaginary (S) parts of P(ix).
    P_ix = P_poly.subs(x, sp.I * x)
    Q_poly = sp.re(P_ix.expand())
    S_poly = sp.im(P_ix.expand())

    # Step 4: Calculate the trace of M1 (Bezout Matrix).
    # M1 = Bm[Q(x), S(x), x]
    # The Bezoutian is defined as (Q(y)S(x) - Q(x)S(y)) / (y-x)
    bezoutian_form = (Q_poly.subs(x, y) * S_poly.subs(x, x) - Q_poly.subs(x, x) * S_poly.subs(x, y)) / (y - x)
    bezoutian_form_expanded = sp.expand(bezoutian_form)

    n_bezout = max(sp.degree(Q_poly, gen=x), sp.degree(S_poly, gen=x))
    
    trace_M1 = 0
    # The trace is the sum of diagonal elements of the Bezout matrix.
    # A diagonal element B[k,k] is the coefficient of x^k * y^k in the Bezoutian.
    for i in range(n_bezout + 1):
        # Extract coefficient of y^i
        coeff_y = bezoutian_form_expanded.coeff(y, i)
        # Extract coefficient of x^i from the result
        coeff_xy = coeff_y.coeff(x, i)
        trace_M1 += coeff_xy

    # Step 5: Calculate the trace of M2 (Sylvester Matrix).
    # M2 = Sm[Q(x), x^10 + S(x), x]
    p1 = Q_poly
    p2 = x**10 + S_poly
    
    # Use the standard formula for the trace of a Sylvester matrix:
    # Tr(Sm(p1, p2)) = n*lc(p1) + m*lc(p2), where m=deg(p1), n=deg(p2), lc is leading coeff.
    deg_p1 = sp.degree(p1, gen=x)
    deg_p2 = sp.degree(p2, gen=x)
    lc_p1 = p1.coeff(x, deg_p1)
    lc_p2 = p2.coeff(x, deg_p2)
    trace_M2 = deg_p2 * lc_p1 + deg_p1 * lc_p2

    # Step 6: Compute the final trace T.
    # T = Tr(M1 (x) I2 + M2) = Tr(M1)*Tr(I2) + Tr(M2)
    trace_I2 = 2
    T = trace_M1 * trace_I2 + trace_M2
    
    print(f"The values of the group properties are R1={R1}, R2={R2}, R3={R3}, R4={R4}.")
    print(f"The constant C is floor(({R1**2}+{R2**2}+{R3**2}+{R4**2})/({R1}+{R2}+{R3}+{R4})) = floor({sum_sq_R}/{sum_R}) = {C}.")
    print(f"The polynomial P(x) is 8*x**9 + 8*x**8 + 8*x**7 + 8*x**6 + 8*x**5 + 8*x**4 + 8*x**3 + 8*x**2 + 8*x.")
    print(f"The real part Q(x) is {Q_poly}.")
    print(f"The imaginary part S(x) is {S_poly}.")
    print(f"The trace of M1 (Bezout Matrix) is {trace_M1}.")
    print(f"The trace of M2 (Sylvester Matrix) is {trace_M2}.")
    print(f"The final trace T is 2*Tr(M1) + Tr(M2) = 2*({trace_M1}) + {trace_M2} = {T}.")

solve()