import sympy as sp

def solve_pendulum_period():
    """
    This function symbolically derives the period of oscillation for the described system
    using Lagrangian mechanics.
    """
    # Step 1: Define symbols
    # R: radius of the disks
    # g: acceleration due to gravity
    # M: mass of each disk
    # t: time
    R, g, M, t = sp.symbols('R g M t', positive=True)

    # The rod has length L = 4R
    L = 4 * R

    # Step 2: Define generalized coordinates and their derivatives
    # x(t): horizontal position of the center of the top disk
    # theta(t): angle of the rod with the vertical
    x = sp.Function('x')(t)
    theta = sp.Function('theta')(t)

    x_dot = sp.diff(x, t)
    theta_dot = sp.diff(theta, t)
    x_ddot = sp.diff(x_dot, t)
    theta_ddot = sp.diff(theta_dot, t)

    # Step 3: Formulate Kinetic and Potential Energy with small angle approximations
    # For small oscillations, we can approximate the full Lagrangian.
    # The resulting Lagrangian in quadratic (in velocities and position) form is:
    # T_approx = 3/2*M*x_dot**2 + M*L*x_dot*theta_dot + 1/2*M*L**2*theta_dot**2
    # V_approx = 1/2*M*g*L*theta**2
    
    # We will use this simplified Lagrangian to derive the equations of motion.
    # We define the equations directly from the derived terms.
    
    print("Deriving the equations of motion for the system...")
    
    # Equation of motion for the 'x' coordinate
    # d/dt(dL/dx_dot) - dL/dx = 0
    # From the Lagrangian, this gives: 3*M*x_ddot + M*L*theta_ddot = 0
    eq_x = 3 * M * x_ddot + M * L * theta_ddot
    print(f"Equation for x: {sp.simplify(eq_x / M)} = 0")

    # Equation of motion for the 'theta' coordinate
    # d/dt(dL/d(theta_dot)) - dL/d(theta) = 0
    # From the Lagrangian, this gives: M*L*x_ddot + M*L**2*theta_ddot + M*g*L*theta = 0
    eq_theta = M * L * x_ddot + M * L**2 * theta_ddot + M * g * L * theta
    print(f"Equation for theta: {sp.simplify(eq_theta / (M*L))} = 0")
    
    # Step 4: Solve the system of equations
    # From the first equation, we can express x_ddot in terms of theta_ddot
    x_ddot_sol = sp.solve(eq_x, x_ddot)[0]
    
    # Substitute this into the second equation to get an equation only in theta
    final_eq = eq_theta.subs(x_ddot, x_ddot_sol)
    
    # The equation is of the form: (Term A) * theta_ddot + (Term B) * theta = 0
    # The angular frequency squared is omega^2 = (Term B) / (Term A)
    term_A = final_eq.coeff(theta_ddot)
    term_B = final_eq.coeff(theta)
    
    omega_squared = term_B / term_A
    
    print(f"\nSolved for the square of the angular frequency, ω²:")
    print(f"ω² = {omega_squared}")

    # Step 5: Calculate the period T = 2*pi / omega
    pi = sp.pi
    period = 2 * pi / sp.sqrt(omega_squared)
    
    # To format the output nicely, we find the numerical coefficients
    term_inside_sqrt = 1 / omega_squared
    num, den = sp.fraction(term_inside_sqrt)
    
    num_coeff = num.coeff(R)
    den_coeff = den.coeff(g)
    
    print("\nThe final equation for the period of motion is:")
    print(f"T = 2 * π * sqrt( ({num_coeff} * R) / ({den_coeff} * g) )")

solve_pendulum_period()