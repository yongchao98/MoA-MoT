import sympy

def solve_bubble_frequency_correction():
    """
    This script calculates the nonlinear frequency correction for the Rayleigh-Plesset equation
    using the Lagrangian/energy method.
    """
    # Define symbolic variables for displacement (x) and polytropic index (gamma)
    x, gamma = sympy.symbols('x gamma')

    # 1. Define Kinetic and Potential Energy terms
    # The system is conservative. The first integral of the R-P equation is:
    # 1/2 * R^3 * R_dot^2 + V(R) = constant
    # where V(R) = R^3/3 - R^(3-3*gamma)/(3-3*gamma)
    # Let R = 1 + x.
    # The effective mass is m(x) = (1+x)^3
    # The potential is V(x)
    m_x = (1 + x)**3
    V_x = (1 + x)**3 / 3 - (1 + x)**(3 - 3*gamma) / (3 - 3*gamma)

    # 2. Expand m(x) and V(x) in Taylor series around x=0
    # The Lagrangian is L = 1/2 * m(x) * x_dot^2 - V(x)
    # m(x) = m0 * (1 + c1*x + c2*x^2 + ...)  (here m0=1)
    # V(x) = V0 + 1/2*k*x^2 + 1/6*a*x^3 + 1/24*b*x^4 + ...

    # Expand mass term to find c1, c2
    m_series = m_x.series(x, 0, 3).removeO()
    c1 = m_series.coeff(x, 1)
    c2 = m_series.coeff(x, 2)

    # Expand potential term to find k, a, b
    V_series = V_x.series(x, 0, 5).removeO()
    k = V_series.coeff(x, 2) * 2
    a = V_series.coeff(x, 3) * 6
    b = V_series.coeff(x, 4) * 24

    # The linearized frequency squared is omega0^2 = k/m0 = k
    omega0_sq = k

    # 3. Use Landau's formula for the frequency correction coefficient k2
    # For a Lagrangian L = 1/2*m0*(1+c1*x+c2*x^2)*x_dot^2 - (1/2*k*x^2 + 1/6*a*x^3 + 1/24*b*x^4)
    # The frequency correction coefficient k2 in omega/omega0 = 1 + k2*epsilon^2 is:
    # k2 = (1/(2*omega0^2)) * (b/(4*m0) - 5*a^2/(12*m0^2*omega0^2)) is not the full formula.
    # The full formula is:
    k2 = (b / (8 * omega0_sq) - 5 * a**2 / (24 * omega0_sq**2) +
          c1 * a / (4 * omega0_sq) - c1**2 / 8 + c2 / 4)

    # 4. Simplify the expression for k2
    k2_simplified = sympy.simplify(k2)

    # The result is a polynomial in gamma.
    k2_poly = sympy.Poly(k2_simplified, gamma)

    # Extract coefficients of the polynomial in ascending order of power
    coeffs = k2_poly.all_coeffs()
    # Sympy returns coeffs in descending order of power, so we reverse them
    coeffs.reverse()
    
    c0 = coeffs[0]
    c1 = coeffs[1]
    c2 = coeffs[2]

    print("The nonlinear frequency correction is given by omega = omega_0 * (1 + k2 * epsilon^2), where epsilon is the amplitude.")
    print(f"The correction coefficient k2 is a polynomial in gamma: k2 = c0 + c1*gamma + c2*gamma^2")
    print("The numbers in the final equation for k2 are:")
    print(f"c0 (constant term) = {c0}")
    print(f"c1 (coefficient of gamma) = {c1}")
    print(f"c2 (coefficient of gamma^2) = {c2}")
    
    # 5. Identify the 3rd term
    # We interpret the "3rd term" as the coefficient of the third term in the
    # polynomial ordered by ascending powers of gamma (i.e., the coefficient of gamma^2).
    third_term_coeff = c2
    
    print("\nThe 3rd term of the nonlinear correction corresponds to the coefficient of the gamma^2 term.")
    print("The value of this term is:")
    print(f"<<<{third_term_coeff}>>>")

solve_bubble_frequency_correction()