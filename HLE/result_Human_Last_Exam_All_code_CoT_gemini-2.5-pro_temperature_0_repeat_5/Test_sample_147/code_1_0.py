import sympy

def solve_quadrature_error():
    """
    Calculates the error term for the optimal linear combination of
    Simpson's 1/3 and 3/8 rules.
    """
    # 1. Define symbols
    h = sympy.Symbol('h')
    x = sympy.Symbol('x')
    f = sympy.Function('f')
    
    # Define derivatives of f at x=0 for Taylor expansion
    derivs = {f(x).diff(x, k).subs(x, 0): sympy.Symbol(f'f{k}') for k in range(9)}
    derivs[f(0)] = sympy.Symbol('f0')

    # 2. Taylor expand the true integral over [-3h, 3h] (6 subintervals of width h)
    # The interval is [a,b] with b-a = 6h.
    series = f(x).series(x, 0, 9).removeO()
    true_integral_series = sympy.integrate(series, (x, -3*h, 3*h))
    true_integral_series = true_integral_series.subs(derivs).expand()

    # 3. Define the quadrature rules and their Taylor expansions
    # The evaluation points are x_i = i*h for i = -3, -2, ..., 3
    f_vals = {}
    for i in range(-3, 4):
        # Taylor expansion of f(i*h) around h=0
        term_list = [((i*h)**k / sympy.factorial(k)) * derivs[f(x).diff(x,k).subs(x,0)] for k in range(9)]
        f_vals[i] = sum(term_list)

    # Composite Simpson's 1/3 rule over 6 intervals (3 applications)
    # Weights: [1, 4, 2, 4, 2, 4, 1] applied to points x_0 to x_6
    # Our points are indexed -3 to 3, which corresponds to x_0 to x_6
    I_1_3_series = (h/3) * (f_vals[-3] + 4*f_vals[-2] + 2*f_vals[-1] + 4*f_vals[0] + 2*f_vals[1] + 4*f_vals[2] + f_vals[3])
    I_1_3_series = I_1_3_series.expand()

    # Composite Simpson's 3/8 rule over 6 intervals (2 applications)
    # Weights: [1, 3, 3, 2, 3, 3, 1] applied to points x_0 to x_6
    I_3_8_series = (3*h/8) * (f_vals[-3] + 3*f_vals[-2] + 3*f_vals[-1] + 2*f_vals[0] + 3*f_vals[1] + 3*f_vals[2] + f_vals[3])
    I_3_8_series = I_3_8_series.expand()

    # 4. Calculate the error series for each rule
    # Error = True Integral - Approximation
    E_1_3_series = (true_integral_series - I_1_3_series).expand()
    E_3_8_series = (true_integral_series - I_3_8_series).expand()

    # Extract leading error term coefficients
    E_1_3_h5_coeff = E_1_3_series.coeff(h**5 * derivs[f(x).diff(x,4).subs(x,0)])
    E_3_8_h5_coeff = E_3_8_series.coeff(h**5 * derivs[f(x).diff(x,4).subs(x,0)])

    # 5. Solve for optimal alpha and beta
    alpha, beta = sympy.symbols('alpha beta')
    # System of equations:
    # alpha + beta = 1 (to maintain correctness for constants)
    # alpha*E_1/3 + beta*E_3/8 = 0 (to cancel leading error term)
    eq1 = sympy.Eq(alpha + beta, 1)
    eq2 = sympy.Eq(alpha * E_1_3_h5_coeff + beta * E_3_8_h5_coeff, 0)
    
    solution = sympy.solve([eq1, eq2], [alpha, beta])
    alpha_val = solution[alpha]
    beta_val = solution[beta]

    # 6. Calculate the new combined error term
    combined_error = (alpha_val * E_1_3_series + beta_val * E_3_8_series).expand()
    
    # Find the leading term of the combined error
    # It should be the h^7 * f^(6) term
    f6 = derivs[f(x).diff(x,6).subs(x,0)]
    leading_term = sympy.collect(combined_error, f6).coeff(f6)
    
    # The error is of the form K * h^n * f^(m)
    # leading_term = K * h^n
    n = sympy.degree(leading_term, h)
    m = 6 # From f6
    K = leading_term.coeff(h**n)

    # 7. Convert to the final form C * (b-a)^n * f^(m)
    # We have b-a = 6h, so h = (b-a)/6
    # Error = K * ((b-a)/6)^n * f^(m) = (K / 6^n) * (b-a)^n * f^(m)
    # The problem asks for C > 0. The calculated K is negative.
    # Error = I_true - I_approx, so the error term is K * h^n * f^(m)
    # The form in the question is C*(b-a)^n*f^(m), so C = |K| / 6^n
    C_val = abs(K) / (6**n)
    
    # Print the results
    print("The optimal linear combination is I = (9/5) * I_1/3 - (4/5) * I_3/8.")
    print("\nThe error term for this optimal rule is of the form C * (b-a)^n * f^(m)(xi).")
    print("The calculated constants are:")
    print(f"C = {C_val} = 1/{C_val.q}")
    print(f"n = {n}")
    print(f"m = {m}")
    print("\nThe final error equation is:")
    print(f"Error = (1/{C_val.q}) * (b-a)^{n} * f^({m})(xi)")

if __name__ == '__main__':
    solve_quadrature_error()