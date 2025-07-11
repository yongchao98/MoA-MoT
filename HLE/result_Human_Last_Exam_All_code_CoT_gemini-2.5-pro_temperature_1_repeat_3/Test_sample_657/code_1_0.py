import sympy

def solve_controller_parametrization():
    """
    This function determines the set of all proper stabilizing controllers H_2(s)
    for the system H_1(s) = s / (s^2 - 1) using Youla-Kucera parametrization.
    """
    # Initialize the symbolic variable 's'
    s = sympy.Symbol('s')
    K = sympy.Function('K')(s)

    # Step 1: Define the plant's numerator and denominator polynomials
    n_p = s
    d_p = s**2 - 1
    
    # We will create a coprime factorization H_1 = N/D where N and D are stable and proper.
    # To do this, we divide n_p and d_p by a stable polynomial m(s) of degree at least 2.
    # A common choice is m(s) = (s+c)^2 for some c > 0. Let's choose c=2.
    c = 2
    m_p = (s + c)**2
    
    # The stable, proper factors are N(s) = n_p/m_p and D(s) = d_p/m_p.
    N = n_p / m_p
    D = d_p / m_p

    # Step 2: Solve the Bezout identity N*X + D*Y = 1 for stable, proper X and Y.
    # This is equivalent to solving the polynomial equation: n_p*X' + d_p*Y' = m_p
    # where X = X'/d_X and Y = Y'/d_Y are the rational functions we seek.
    # A systematic way to find a proper and stable solution is to solve the Diophantine equation
    # s * nX(s) + (s^2-1) * nY(s) = (s+c)^3 for polynomials nX and nY,
    # where X(s) = nX(s)/(s+c) and Y(s) = nY(s)/(s+c).
    
    # To ensure X(s) is proper, deg(nX) must be <= 1. To ensure Y(s) is proper, deg(nY) must be <= 1.
    # From the equation s*nX + (s^2-1)*nY = (s+c)^3, by comparing degrees, we find that
    # deg(nY) must be 1. Let nY = a*s + b.
    a, b = sympy.symbols('a b')
    nY_poly = a*s + b
    
    # Substitute nY into the equation and solve for nX.
    # s*nX = (s+c)**3 - (s**2-1)*nY
    rhs = (s+c)**3 - (s**2 - 1) * nY_poly
    
    # For s*nX to be a polynomial divisible by s, its constant term must be zero.
    const_term = rhs.subs(s, 0)
    # const_term = c**3 - (-1)*b = c**3 + b = 0 => b = -c**3
    b_val = -c**3
    
    # Substitute b back into the equation for the right-hand side.
    rhs = rhs.subs(b, b_val)
    
    # Now, nX is rhs/s
    nX_poly_gen = sympy.simplify(rhs / s)

    # For X(s) = nX/(s+c) to be proper, deg(nX) must be at most 1.
    # This means the coefficient of s^2 in nX must be zero.
    s_squared_coeff = sympy.Poly(nX_poly_gen, s).coeff_monomial(s**2)
    # s_squared_coeff = 1-a = 0 => a = 1
    a_val = sympy.solve(s_squared_coeff, a)[0]

    # Now we have the final polynomials nX and nY.
    nX = sympy.simplify(nX_poly_gen.subs(a, a_val))
    nY = sympy.simplify(nY_poly.subs({a: a_val, b: b_val}))
    
    # So, the particular solution to the Bezout identity is:
    X = nX / (s + c)
    Y = nY / (s + c)

    # Step 3: Construct the general formula for the controller H_2(s).
    # H_2(s) = (X + D*K) / (Y - N*K)
    # To simplify, multiply numerator and denominator by m_p = (s+c)^2
    # H_2(s) = (X*m_p + D*m_p*K) / (Y*m_p - N*m_p*K)
    # H_2(s) = (X*(s+c)^2 + d_p*K) / (Y*(s+c)^2 - n_p*K)
    # H_2(s) = ( (nX/(s+c))*(s+c)^2 + d_p*K ) / ( (nY/(s+c))*(s+c)^2 - n_p*K )
    # H_2(s) = ( nX*(s+c) + d_p*K ) / ( nY*(s+c) - n_p*K )

    num_poly_part = sympy.expand(nX * (s + c))
    den_poly_part = sympy.expand(nY * (s + c))

    # The final expression for the set of all proper stabilizing controllers H_2(s)
    # parametrized by a stable proper function K(s).
    
    # The prompt requests the output in a specific format, showing each term.
    # We will build the string for the numerator and denominator.
    
    num_str = f"({d_p})*K(s) + ({num_poly_part})"
    den_str = f"({-n_p})*K(s) + ({den_poly_part})"
    
    final_eq = f"H_2(s) = {num_str} / {den_str}"

    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print(final_eq)
    print("\nwhere K(s) is any stable and proper rational function.")
    
    # For clarity, let's print the expanded version with coefficients
    num_exp_poly = sympy.Poly(num_poly_part, s)
    den_exp_poly = sympy.Poly(den_poly_part, s)
    
    c2, c1, c0 = num_exp_poly.all_coeffs()
    d2, d1, d0 = den_exp_poly.all_coeffs()
    
    print("\nExpanded form:")
    print(f"H_2(s) = (({s**2} - 1)*K(s) + {c2}*s**2 + {c1}*s + {c0}) / (-{s}*K(s) + {d2}*s**2 + {d1}*s + {d0})")


solve_controller_parametrization()
<<<H_2(s) = ((s**2 - 1)*K(s) + 14*s**2 + 41*s + 26) / (-s*K(s) + s**2 - 6*s - 16)>>>