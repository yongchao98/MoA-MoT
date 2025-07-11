import sympy as sp

def solve_theta_prime():
    """
    This function symbolically calculates the expression for theta'(t).
    """
    # Define the symbols used in the derivation.
    # t: time
    # c: a given constant
    # K: Gaussian curvature, a function of t
    # theta: the angle, a function of t
    # A, B: coefficients, functions of t
    t = sp.Symbol('t')
    c = sp.Symbol('c', real=True, positive=True)
    K = sp.Function('K')(t)
    theta = sp.Function('theta')(t)
    A = sp.Function('A')(t)
    B = sp.Function('B')(t)

    # From the physical setup, we derived the system of ODEs for A and B:
    # A' = -(K/c) * B
    # B' = c * A
    A_prime = -(K / c) * B
    B_prime = c * A

    # The rate of change of the angle theta is given by the formula:
    # theta' = (A*B' - B*A') / (A^2 + B^2)
    theta_prime_numerator = A * B_prime - B * A_prime
    theta_prime_denominator = A**2 + B**2
    theta_prime_expr = theta_prime_numerator / theta_prime_denominator

    # We now substitute the polar coordinate representations for A and B
    # to express the result in terms of theta.
    # R(t) is the radial part, which will cancel out.
    R = sp.Function('R')(t)
    A_polar = R * sp.cos(theta)
    B_polar = R * sp.sin(theta)

    # Substitute the polar forms of A and B into the expression for theta'.
    theta_prime_final = theta_prime_expr.subs({A: A_polar, B: B_polar})

    # Simplify the expression.
    theta_prime_simplified = sp.simplify(theta_prime_final)

    # To make the output match the format of the answers, we replace
    # the function forms K(t) and theta(t) with simple symbols.
    K_sym = sp.Symbol('K(gamma(t))')
    theta_sym = sp.Symbol('theta(t)')
    cos_sym = sp.Symbol('cos(theta(t))')
    sin_sym = sp.Symbol('sin(theta(t))')

    final_readable_expr = theta_prime_simplified.subs({
        K: K_sym,
        sp.cos(theta): cos_sym,
        sp.sin(theta): sin_sym
    })
    
    # We add the numbers and operators explicitly as requested for clarity.
    cos_squared = sp.Pow(cos_sym, 2)
    sin_squared = sp.Pow(sin_sym, 2)

    term1 = c * cos_squared
    term2 = (sp.S(1)/c) * K_sym * sin_squared

    final_equation_form = sp.Add(term1, term2, evaluate=False)

    print("The final expression for theta'(t) is:")
    print(f"{c} * {cos_squared} + ({sp.S(1)}/{c}) * {K_sym} * {sin_squared}")


solve_theta_prime()