import sympy

def solve_pendulum_period():
    """
    This function calculates and displays the period of the pendulum-like system
    using symbolic mathematics.
    """
    # Define the symbolic variables for the problem
    R, g, t = sympy.symbols('R g t', positive=True)
    # Define theta as a function of time for the equation of motion
    theta = sympy.Function('theta')(t)

    # From the Lagrangian derivation, we obtain the linearized equation of motion for small theta.
    # The equation has the form: (effective mass term) * d^2(theta)/dt^2 + (restoring force term) * theta = 0

    # The coefficient for the second derivative of theta is 8*R/3
    mass_coeff = sympy.S(8) * R / 3
    # The coefficient for theta (from the potential energy part) is g
    spring_coeff = g

    # The full equation of motion is:
    eom = sympy.Eq(mass_coeff * theta.diff(t, t) + spring_coeff * theta, 0)
    
    print("The derived linearized equation of motion for small oscillations is:")
    sympy.pprint(eom)
    print("\nThis equation is in the standard form of a simple harmonic oscillator: M_eff * d^2(theta)/dt^2 + K_eff * theta = 0")
    print("Where the effective 'mass' M_eff is the coefficient of the d^2(theta)/dt^2 term, and the effective 'spring constant' K_eff is the coefficient of the theta term.")
    
    # We can rewrite it in the more common form: d^2(theta)/dt^2 + omega^2 * theta = 0
    # where omega^2 = K_eff / M_eff
    omega_squared = spring_coeff / mass_coeff

    print(f"\nFrom this, the square of the angular frequency (omega^2) is K_eff / M_eff:")
    omega_sq_eq = sympy.Eq(sympy.Symbol('omega')**2, omega_squared)
    sympy.pprint(omega_sq_eq)

    # The period of motion is T = 2*pi / omega
    period = 2 * sympy.pi / sympy.sqrt(omega_squared)

    print("\nThe period of motion (P) is 2*pi / omega, which is 2*pi / sqrt(omega^2):")
    period_eq = sympy.Eq(sympy.Symbol('P'), period)
    sympy.pprint(period_eq)

solve_pendulum_period()