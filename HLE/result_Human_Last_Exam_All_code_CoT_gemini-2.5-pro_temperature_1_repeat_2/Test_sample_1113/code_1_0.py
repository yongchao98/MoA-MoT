import math

def solve_robot_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters
    h = 1.0  # Robot height in meters
    r = 0.25  # Arm length in meters
    l_c = 10.0  # Chain length in meters
    v = 10.0  # Robot speed in m/s
    dot_beta = 1.0  # Arm angular speed in rad/s
    d = 20.0  # Visible path length (diameter) in meters
    l_shadow = 10.0 * math.sqrt(3)  # Shadow length in meters

    # Step 1: Calculate geometric parameters
    R = d / 2.0  # Radius of the circular path
    # The tilt angle theta of the path plane with the ground
    cos_theta = l_shadow / d
    theta = math.acos(cos_theta)
    sin_theta = math.sin(theta)

    # Step 2: Determine motion parameters
    # The robot's angular speed along the path
    omega_robot = v / R
    
    # Based on the problem statement, omega_robot = 10/10 = 1 rad/s.
    # The arm's rotation angle is beta(t) = dot_beta * t = 1 * t = t.
    # The robot's angular position is alpha(t) = alpha_P + omega_robot * t = alpha_P + t.

    # Step 3: Set the starting position alpha_P
    # The problem states the highest point (alpha=pi) is reached in the first quarter of the journey (t <= T/4 = pi/2).
    # This implies pi/2 <= alpha_P <= pi.
    # We assume the boundary case where the journey starts at alpha_P = pi/2.
    alpha_P = math.pi / 2.0

    # Step 4: Formulate and solve the equation for t
    # The height of the arm's tip is z_tip(t). The chain loses contact when z_tip(t) = l_c.
    # z_tip(t) = R*(1-cos(alpha_P+t))*sin_theta + h*cos_theta + r*cos(t)*sin(alpha_P+t)*sin_theta + r*sin(t)*cos_theta
    # With alpha_P = pi/2, cos(pi/2+t) = -sin(t) and sin(pi/2+t) = cos(t).
    # l_c = R*(1+sin(t))*sin_theta + h*cos_theta + r*cos(t)*cos(t)*sin_theta + r*sin(t)*cos_theta
    # Rearranging this gives a quadratic equation in sin(t): a*x^2 + b*x + c = 0 where x = sin(t).

    # Coefficients of the quadratic equation 0.125*x^2 - (5 + sqrt(3)/8)*x + (4.875 - sqrt(3)/2) = 0
    # Let's write it as a*x^2 + b*x + c = 0
    a = -r * sin_theta 
    b = R * sin_theta + r * cos_theta
    c = R * sin_theta + h * cos_theta - l_c
    
    # Substituting cos^2(t) = 1 - sin^2(t), the equation becomes:
    # l_c = R*(1+sin(t))*sin_theta + h*cos_theta + r*(1-sin^2(t))*sin_theta + r*sin(t)*cos_theta
    # l_c = R*sin_theta + R*sin_theta*sin(t) + h*cos_theta + r*sin_theta - r*sin_theta*sin^2(t) + r*cos_theta*sin(t)
    # (-r*sin_theta)*sin^2(t) + (R*sin_theta + r*cos_theta)*sin(t) + (R*sin_theta + h*cos_theta + r*sin_theta - l_c) = 0
    
    quad_a = -r * sin_theta
    quad_b = R * sin_theta + r * cos_theta
    quad_c = R * sin_theta + h * cos_theta + r * sin_theta - l_c
    
    print("The problem reduces to solving a quadratic equation for sin(t) of the form: a*(sin(t))^2 + b*sin(t) + c = 0")
    print(f"The calculated coefficients are:")
    print(f"a = {quad_a:.4f}")
    print(f"b = {quad_b:.4f}")
    print(f"c = {quad_c:.4f}")
    
    # Solve the quadratic equation for sin_t
    discriminant = quad_b**2 - 4 * quad_a * quad_c
    
    if discriminant < 0:
        print("No real solution for time exists.")
        return

    # We expect two solutions for sin_t
    sin_t1 = (-quad_b + math.sqrt(discriminant)) / (2 * quad_a)
    sin_t2 = (-quad_b - math.sqrt(discriminant)) / (2 * quad_a)

    # Choose the physically valid solution for sin_t (must be between -1 and 1)
    valid_sin_t = None
    if -1 <= sin_t1 <= 1:
        valid_sin_t = sin_t1
    if -1 <= sin_t2 <= 1:
        # We need the smallest positive time t, which corresponds to the smaller positive sin(t)
        if valid_sin_t is None or sin_t2 < valid_sin_t:
             valid_sin_t = sin_t2

    if valid_sin_t is None:
        print("No valid solution for sin(t) found.")
        return

    # Calculate time t (smallest positive value)
    t = math.asin(valid_sin_t)

    print("\nThe final equation to solve is:")
    print(f"{quad_a:.4f} * (sin(t))^2 + {quad_b:.4f} * sin(t) + {quad_c:.4f} = 0")
    print("\nSolving for the smallest positive t:")
    print(f"The chain will first lose contact with the ground at t = {t:.2f} seconds.")


solve_robot_chain_problem()