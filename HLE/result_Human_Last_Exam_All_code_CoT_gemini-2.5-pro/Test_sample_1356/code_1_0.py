import sympy
from sympy import symbols, Function, sin, cos, diff, Eq, solve, simplify, pi, sqrt

def solve_pendulum_period():
    """
    This function symbolically derives the period of oscillation for the described system
    using Lagrangian mechanics.
    """
    print("Deriving the period of motion step-by-step using Lagrangian mechanics.")

    # 1. Define symbols and generalized coordinates
    M, R, g, t = symbols('M R g t', positive=True)
    theta = Function('theta')(t)
    x = Function('x')(t)

    # Define first and second derivatives
    theta_dot = diff(theta, t)
    theta_ddot = diff(theta_dot, t)
    x_dot = diff(x, t)
    x_ddot = diff(x_dot, t)

    print("\nStep 1: Defining Kinetic and Potential Energy")
    # 2. Calculate Kinetic Energy (T)
    # KE of Disk 1 (top, rolling without slipping)
    # I_disk = 1/2 * M * R**2
    # omega_1 = x_dot / R
    # T1_trans = 1/2 * M * x_dot**2
    # T1_rot = 1/2 * I_disk * omega_1**2 = 1/2 * (1/2 * M * R**2) * (x_dot/R)**2 = 1/4 * M * x_dot**2
    T1 = sympy.Rational(3, 4) * M * x_dot**2

    # KE of Disk 2 (hanging)
    # Position of center of Disk 2: (x + 4*R*sin(theta), -4*R*cos(theta))
    # Velocity of center of Disk 2
    x2_dot = diff(x + 4 * R * sin(theta), t)
    y2_dot = diff(-4 * R * cos(theta), t)
    # T2_trans = 1/2 * M * (x2_dot**2 + y2_dot**2)
    T2_trans = simplify(sympy.Rational(1, 2) * M * (x2_dot**2 + y2_dot**2))
    
    # Rotational KE of Disk 2 (rotates with the rod, angular velocity theta_dot)
    T2_rot = sympy.Rational(1, 2) * (sympy.Rational(1, 2) * M * R**2) * theta_dot**2
    
    T = T1 + T2_trans + T2_rot
    T = simplify(T)

    # 3. Calculate Potential Energy (V)
    # V1 = 0 (on the table)
    # V2 = M * g * y2 = M * g * (-4*R*cos(theta))
    V = -4 * M * g * R * cos(theta)

    print(f"Total Kinetic Energy (T): {T}")
    print(f"Total Potential Energy (V): {V}")
    
    # 4. Formulate the Lagrangian
    L = T - V

    print("\nStep 2: Deriving Equations of Motion from the Lagrangian")
    # 5. Derive Lagrange's equations
    # Equation for x
    # Since dL/dx = 0, the conjugate momentum is conserved.
    # We assume the system starts from rest, so this constant is 0.
    px = diff(L, x_dot)
    # px = 0 for a system starting from rest
    constraint_eq = Eq(px, 0)
    
    # Solve for x_dot in terms of theta_dot for small angles (cos(theta) -> 1)
    x_dot_sol = solve(constraint_eq.subs(cos(theta), 1), x_dot)[0]
    
    print(f"From conservation of momentum in x: {px} = 0")
    print(f"For small angles, this gives a constraint: x_dot = {x_dot_sol}")
    
    # Equation for theta
    L_theta_eq = simplify(diff(diff(L, theta_dot), t) - diff(L, theta))

    print("\nStep 3: Linearizing for Small Oscillations")
    # 6. Linearize the equation for theta
    # Substitute the constraint for x_ddot
    x_ddot_sol = diff(x_dot_sol, t)
    L_theta_eq_sub = L_theta_eq.subs(x_ddot, x_ddot_sol)

    # Apply small angle approximations: sin(theta) -> theta, cos(theta) -> 1
    # and drop higher order terms (like theta_dot**2, which arise from differentiation)
    # Manually linearizing by substituting approximations into the equation
    # A more robust way is to expand L as a series, but this is sufficient.
    linear_eq = L_theta_eq_sub.subs([(sin(theta), theta), (cos(theta), 1), (theta_dot**2, 0)]).doit()
    linear_eq = simplify(linear_eq)
    
    final_eq = Eq(linear_eq, 0)
    print(f"The linearized equation of motion for theta is: {final_eq}")
    
    # 7. Solve for angular frequency squared (omega^2)
    # The equation is of the form: I_eff * theta_ddot + k_eff * theta = 0
    # where omega^2 = k_eff / I_eff
    rearranged_eq = solve(final_eq, theta_ddot)[0]
    omega_sq = -rearranged_eq / theta
    
    print("\nStep 4: Calculating the Period")
    print(f"The equation is in the form of a simple harmonic oscillator: theta_ddot + omega^2 * theta = 0")
    print(f"Angular frequency squared (omega^2) = {omega_sq}")

    # 8. Calculate the period T
    period = 2 * pi * sqrt(1 / omega_sq)
    
    # Extracting numerical coefficients for clear output
    num, den = omega_sq.as_numer_denom()
    num_val = num / (g)
    den_val = den / (R)

    print("\nThe final expression for the period of motion is:")
    print(f"T = 2 * pi * sqrt({den_val}*R / ({num_val}*g))")

solve_pendulum_period()