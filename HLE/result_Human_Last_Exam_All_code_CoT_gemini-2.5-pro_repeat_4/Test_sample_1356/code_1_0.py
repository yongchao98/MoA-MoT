import sympy
from sympy import symbols, Function, Eq, sin, cos, diff, solve, sqrt, pi

def solve_pendulum_period():
    """
    This function uses symbolic mathematics to derive the period of oscillation
    for the described two-disk pendulum system.
    """
    # Step 1: Define symbols and functions for our generalized coordinates
    m, R, g, t = symbols('m R g t', positive=True)
    x = Function('x')(t)      # Position of the center of the top disk
    theta = Function('theta')(t)  # Angle of the rod with the vertical

    # Define time derivatives (dots)
    x_dot = diff(x, t)
    x_ddot = diff(x_dot, t)
    theta_dot = diff(theta, t)
    theta_ddot = diff(theta_dot, t)

    # Step 2: Define Kinetic and Potential Energy
    
    # --- Kinetic Energy (T) ---
    # Disk 1 (top disk)
    # Moment of inertia of a solid disk
    I1 = 0.5 * m * R**2
    # Rolling without slipping condition: angular velocity omega1 = x_dot / R
    omega1 = x_dot / R
    # KE of Disk 1 = Translational KE + Rotational KE
    T1 = 0.5 * m * x_dot**2 + 0.5 * I1 * omega1**2
    # Simplify T1
    T1 = sympy.simplify(T1) # Should be (3/4)*m*x_dot**2

    # Disk 2 (bottom disk)
    # Position of Disk 2's center
    x2 = x + 4 * R * sin(theta)
    y2 = R - 4 * R * cos(theta)
    # Velocity of Disk 2's center
    x2_dot = diff(x2, t)
    y2_dot = diff(y2, t)
    # KE of Disk 2 = Translational KE. It does not rotate on its own axis.
    T2 = 0.5 * m * (x2_dot**2 + y2_dot**2)
    # Simplify T2
    T2 = sympy.simplify(T2)

    # Total Kinetic Energy
    T = T1 + T2

    # --- Potential Energy (U) ---
    # U for Disk 1 (height R from a y=0 table)
    U1 = m * g * R
    # U for Disk 2 (height y2)
    U2 = m * g * y2
    # Total Potential Energy
    U = U1 + U2

    # --- Lagrangian (L = T - U) ---
    L = T - U

    # Step 3: Linearize for small angles (theta)
    # For small theta: sin(theta) -> theta, cos(theta) -> 1 - theta^2/2
    # We apply this to the Lagrangian. We also ignore terms of order > 2.
    L_lin = L.subs([(sin(theta), theta), (cos(theta), 1 - theta**2/2)]).doit()
    # Keep terms up to second order in theta and its derivative
    # This can be done by expanding around theta=0, theta_dot=0
    # For this system, manual substitution is easier.
    
    # We will derive the full equations and then linearize them.
    # Lagrange's Equations: d/dt(dL/d(q_dot)) - dL/dq = 0

    # Equation for coordinate x
    eq_x = diff(diff(L, x_dot), t) - diff(L, x)
    eq_x_lin = sympy.limit(eq_x / (m*R), theta, 0).doit()
    
    # Equation for coordinate theta
    eq_theta = diff(diff(L, theta_dot), t) - diff(L, theta)
    # Linearizing requires careful handling of sin(theta)/theta etc. for small theta
    eq_theta_lin = eq_theta.subs([(sin(theta), theta), (cos(theta), 1), (theta_dot**2, 0)]).doit()
    eq_theta_lin = eq_theta_lin / (m * R)

    # Step 4: Solve the system of linear differential equations
    # We have two equations for x_ddot and theta_ddot
    # eq_x_lin = 5/2 * x_ddot/R + 4 * theta_ddot = 0
    # eq_theta_lin = 4 * x_ddot + 16*R*theta_ddot + 4*g*theta = 0
    
    # From the first equation, express x_ddot in terms of theta_ddot
    sol_x_ddot = solve(eq_x_lin, x_ddot)[0]

    # Substitute this into the second equation
    final_eq = eq_theta_lin.subs(x_ddot, sol_x_ddot).simplify()
    # The result is an equation of the form: A*theta_ddot + B*theta = 0
    # This is SHM: theta_ddot + omega^2 * theta = 0
    # So, omega^2 = - (coefficient of theta) / (coefficient of theta_ddot)
    
    # Let's isolate the coefficients
    # final_eq is of the form: 12*R*theta_ddot(t)/5 + g*theta(t) = 0
    p = sympy.Poly(final_eq.lhs, theta_ddot)
    coeff_theta_ddot = p.coeff_monomial(theta_ddot)
    
    # To get the theta coefficient, we set theta_ddot to 0 in the poly
    coeff_theta = p.subs(theta_ddot, 0) / theta

    # Step 5: Calculate angular frequency squared (omega^2)
    omega_squared = -coeff_theta / coeff_theta_ddot
    
    # Step 6: Calculate the Period T
    period = 2 * pi * sqrt(1 / omega_squared)
    
    # --- Output the results ---
    print("The derivation leads to the equation for simple harmonic motion:")
    print(f"{final_eq.lhs} = 0")
    print("\nFrom this, the angular frequency squared (ω²) is found to be:")
    sympy.pprint(Eq(symbols('omega')**2, omega_squared))
    
    print("\nThe period of motion T = 2π/ω is therefore:")
    # Get the numerator and denominator to print the numbers clearly
    num, den = sympy.fraction(1/omega_squared)
    num_val = num.args[0]
    den_val = den.args[0]
    
    # We want to print T = 2*pi*sqrt(12*R/(5*g))
    # 1/omega_squared = 12*R/(5*g)
    num_T_sq = sympy.Poly(num, R).coeffs()[0] # should be 12
    den_T_sq_g = sympy.Poly(den, g).coeffs()[0] # should be 5

    print(f"T = 2 * pi * sqrt({num_T_sq} * R / ({den_T_sq_g} * g))")

solve_pendulum_period()