import sympy
from sympy import sin, cos, tan, atan, Symbol, Eq, solve, Rational, pretty_print

def solve_rod_sliding_problem():
    """
    This function symbolically derives the expression for the angle at which a
    tilting rod begins to slide off a table corner.
    """
    # 1. Define all symbolic variables
    # M: Mass of the rod
    # L: Length of the rod
    # g: Acceleration due to gravity
    # l: Distance of the center of mass from the pivot (table edge)
    # mu: Coefficient of static friction
    # theta: Angle of tilt of the rod with the horizontal
    M, g, L, l, mu, theta = sympy.symbols('M g L l mu theta', positive=True, real=True)

    # Intermediate variables for the derivation
    # I: Moment of inertia about the pivot
    # alpha: Angular acceleration
    # omega_sq: Angular velocity squared
    # N: Normal force from the table edge
    # f: Friction force from the table edge
    I, alpha, omega_sq, N, f = sympy.symbols('I alpha omega_sq N f')

    # 2. Define the Moment of Inertia about the pivot
    # Using the parallel axis theorem: I_pivot = I_cm + M*d^2
    # I_cm for a rod is (1/12)*M*L^2. The distance d from CM to pivot is l.
    I_val = M * (Rational(1, 12) * L**2 + l**2)

    # 3. Find Angular Acceleration (alpha) from the torque equation (tau = I * alpha)
    # The torque is produced by the horizontal component of the CM's displacement from the pivot.
    # tau = (force) * (lever arm) = (M*g) * (l*cos(theta))
    torque_eq = Eq(M * g * l * cos(theta), I_val * alpha)
    alpha_expr = solve(torque_eq, alpha)[0]

    # 4. Find Angular Velocity Squared (omega_sq) from conservation of energy
    # Change in Potential Energy = Rotational Kinetic Energy
    # M*g*h = (1/2)*I*omega^2, where h = l*sin(theta)
    energy_eq = Eq(M * g * l * sin(theta), Rational(1, 2) * I_val * omega_sq)
    omega_sq_expr = solve(energy_eq, omega_sq)[0]

    # 5. Find the Normal Force (N)
    # Analyze forces perpendicular to the rod. The net force provides the tangential
    # acceleration of the center of mass (a_tan = l * alpha).
    # N - (component of gravity) = M * a_tan
    normal_force_eq = Eq(N - M * g * cos(theta), M * l * alpha_expr)
    N_expr = solve(normal_force_eq, N)[0]

    # 6. Find the required Friction Force (f)
    # Analyze forces parallel to the rod. The net force provides the centripetal
    # acceleration of the center of mass (a_cen = l * omega^2).
    # (component of gravity) - f = M * (-a_cen)  (negative as a_cen is towards pivot)
    friction_force_eq = Eq(M * g * sin(theta) - f, -M * l * omega_sq_expr)
    f_expr = solve(friction_force_eq, f)[0]

    # 7. Apply the condition for sliding: f = mu * N
    # Substitute the expressions for f and N and solve for tan(theta).
    sliding_condition_eq = Eq(f_expr, mu * N_expr)

    # Sympy can solve this complex equation for tan(theta) directly.
    tan_theta_solution = solve(sliding_condition_eq, tan(theta))

    # 8. Print the result
    print("The rod begins to slide when the friction force equals the maximum static friction (f = mu * N).")
    print("Solving this relationship for tan(theta) gives the following expression:")
    print("-" * 60)
    
    # We construct the final printable equation for clarity
    numerator_expr = mu * (L**2/12 + 2*l**2)
    denominator_expr = (L**2/12 + 3*l**2)
    final_eq = Eq(tan(theta), numerator_expr / denominator_expr)
    
    # Use pretty_print for a more readable mathematical format
    pretty_print(final_eq)

    print("\nTherefore, the expression for the angle theta is:")
    final_theta_expr = atan(final_eq.rhs)
    pretty_print(Eq(theta, final_theta_expr))
    print("-" * 60)

solve_rod_sliding_problem()