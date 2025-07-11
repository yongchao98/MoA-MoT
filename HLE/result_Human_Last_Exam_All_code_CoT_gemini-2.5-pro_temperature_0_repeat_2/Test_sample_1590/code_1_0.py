import sympy as sp

def solve_rod_sliding_angle():
    """
    This script symbolically derives the expression for the angle theta at which a
    tilting rod begins to slide off the corner of a table.
    """
    # Step 1: Define all symbolic variables
    # L: total length of the rod
    # l: distance from the pivot to the center of mass
    # M: mass of the rod
    # g: acceleration due to gravity
    # mu: coefficient of static friction
    # theta: angle of the rod with the horizontal
    L, l, M, g, mu, theta = sp.symbols('L l M g mu theta', positive=True)

    # Step 2: Define moments of inertia
    # I_cm: Moment of inertia about the center of mass for a thin rod
    I_cm = sp.Rational(1, 12) * M * L**2
    # I_p: Moment of inertia about the pivot point 'P' using the parallel axis theorem
    # The distance from the CM to the pivot is 'l'.
    I_p = I_cm + M * l**2

    # Step 3: Use energy and torque to find angular velocity and acceleration
    # From conservation of energy (Loss in PE = Gain in KE): M*g*l*sin(theta) = 0.5 * I_p * thetadot**2
    # This gives thetadot**2 (angular velocity squared)
    thetadot_sq = (2 * M * g * l * sp.sin(theta)) / I_p

    # From the torque equation (Torque = I_p * alpha): Torque = (M*g*sin(theta)) * l
    # This gives thetaddot (angular acceleration)
    thetaddot = (M * g * l * sp.sin(theta)) / I_p

    # Step 4: Set up the force equations in the rotating frame (parallel and perpendicular to the rod)
    # Perpendicular to the rod: N - M*g*cos(theta) = M * a_tangential = M * l * thetaddot
    N = M * g * sp.cos(theta) + M * l * thetaddot

    # Parallel to the rod: M*g*sin(theta) - f = M * a_radial = M * (-l * thetadot**2)
    # (The radial acceleration is centripetal, pointing towards the pivot, hence negative)
    f = M * g * sp.sin(theta) + M * l * thetadot_sq

    # Step 5: Apply the condition for sliding: f = mu * N
    sliding_condition = sp.Eq(f, mu * N)

    # Step 6: Solve the equation for tan(theta)
    # We can simplify the equation by dividing all terms by M*g and rearranging
    # to solve for tan(theta). Sympy can do this for us.
    # We create a symbol for tan(theta) to make the solver's job easier.
    tan_theta = sp.Symbol('tan(theta)')
    
    # Substitute cos(theta) = 1/sqrt(1+tan(theta)**2) and sin(theta) = tan(theta)/sqrt(1+tan(theta)**2)
    # A simpler way is to rearrange the equation manually and then solve.
    # f = mu * N
    # M*g*sin(theta) + M*l*thetadot_sq = mu * (M*g*cos(theta) + M*l*thetaddot)
    # Substitute expressions for thetadot_sq and thetaddot
    # M*g*sin(theta) + M*l*(2*M*g*l*sin(theta)/I_p) = mu * (M*g*cos(theta) + M*l*(M*g*l*sin(theta)/I_p))
    # Divide by M*g
    # sin(theta) + 2*M*l**2*sin(theta)/I_p = mu * (cos(theta) + M*l**2*sin(theta)/I_p)
    # Divide by cos(theta)
    # tan(theta) + 2*M*l**2*tan(theta)/I_p = mu * (1 + M*l**2*tan(theta)/I_p)
    # tan(theta) * (1 + 2*M*l**2/I_p - mu*M*l**2/I_p) = mu
    # tan(theta) = mu / (1 + (M*l**2*(2-mu))/I_p)
    # tan(theta) = mu * I_p / (I_p + M*l**2*(2-mu))
    
    tan_theta_expr = (mu * I_p) / (I_p + M * l**2 * (2 - mu))
    
    # Substitute I_p back to get the full expression in terms of L and l
    tan_theta_final = tan_theta_expr.subs(I_p, M * (L**2 / 12 + l**2))
    
    # Simplify the expression
    tan_theta_final_simplified = sp.simplify(tan_theta_final)

    # Print the final result
    print("The expression for the tangent of the angle at which the rod begins to slide is:")
    # Using sp.pretty_print for a more readable equation format
    sp.pretty_print(sp.Eq(sp.tan(theta), tan_theta_final_simplified))
    
    print("\nTherefore, the angle theta is:")
    final_expression = sp.atan(tan_theta_final_simplified)
    sp.pretty_print(sp.Eq(theta, final_expression))

solve_rod_sliding_angle()