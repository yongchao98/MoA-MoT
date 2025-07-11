import math

def solve_robot_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    beta_dot = 1.0  # rad/s, arm angular speed
    d = 20.0  # m, visible diameter of the path
    l_shadow = 10.0 * math.sqrt(3)  # m, shadow length

    # 1. Calculate path parameters
    R = d / 2.0  # Radius of the circular path
    # The tilt angle alpha is found from the projection
    # 2*R*cos(alpha) = l_shadow => cos(alpha) = l_shadow / d
    alpha = math.acos(l_shadow / d)  # Tilt angle of the path
    # The center of the circle's height z_c is R*sin(alpha)
    # so the lowest point at z = z_c - R*sin(alpha) is 0.
    z_c = R * math.sin(alpha)

    # Robot's angular speed on the path
    theta_dot = v / R
    # Since v=10 and R=10, theta_dot = 1 rad/s. So theta(t) = t.
    # Arm angle beta(t) = beta_0 + beta_dot*t.
    # beta_0 is pi/2 (forward), beta_dot = 1. So beta(t) = pi/2 + t.

    # 2. Formulate the equation for the z-coordinate of the arm tip
    # We assume the arm is attached to the base of the robot (P_robot)
    # as it's the only interpretation that yields a solution.
    # P_tip_z(t) = z_robot(t) + z_arm(t)
    # z_robot(t) = z_c + R * sin(theta(t)) * sin(alpha)
    # z_arm(t) = r * (cos(beta(t))*n_z + sin(beta(t))*u_forward_z(t))
    # n_z = -cos(alpha)
    # u_forward_z(t) = cos(theta(t))*sin(alpha)
    # Substituting theta=t and beta=pi/2+t:
    # z_robot(t) = 5 + 5*sin(t)
    # z_arm(t) = r * (-sin(t)*(-cos(alpha)) + cos(t)*(cos(t)*sin(alpha)))
    # z_arm(t) = r * (sin(t)*cos(alpha) + cos^2(t)*sin(alpha))
    # P_tip_z(t) = 5 + 5*sin(t) + r*(sin(t)*cos(alpha) + cos^2(t)*sin(alpha))

    # 3. Set P_tip_z(t) = l_c and solve for s = sin(t)
    # 5 + 5s + r*s*cos(alpha) + r*(1-s^2)*sin(alpha) = 10
    # -5 + (5 + r*cos(alpha))*s + r*sin(alpha) - r*sin(alpha)*s^2 = 0
    # (r*sin(alpha))*s^2 - (5 + r*cos(alpha))*s + (5 - r*sin(alpha)) = 0

    # Coefficients of the quadratic equation a*s^2 + b*s + c = 0
    a = r * math.sin(alpha)
    b = -(5.0 + r * math.cos(alpha))
    c = 5.0 - r * math.sin(alpha)
    
    # The problem asks to output the numbers in the final equation.
    # Let's use the simplified equation s^2 - (40 + sqrt(3))s + 39 = 0
    # which is derived from the same physics but with fractions.
    # a_s = 1, b_s = -(40+sqrt(3)), c_s = 39
    a_s = 1.0
    b_s = -(40.0 + math.sqrt(3))
    c_s = 39.0

    print("The problem is solved by finding the root of a quadratic equation for s = sin(t):")
    print(f"a*s^2 + b*s + c = 0")
    print(f"where a = {a_s:.4f}")
    print(f"where b = {b_s:.4f}")
    print(f"where c = {c_s:.4f}")
    
    # Solve the quadratic equation for s = sin(t)
    discriminant = b_s**2 - 4 * a_s * c_s
    
    if discriminant < 0:
        print("The chain never loses contact with the ground.")
        solution_t = None
    else:
        # Two potential solutions for s
        s1 = (-b_s + math.sqrt(discriminant)) / (2 * a_s)
        s2 = (-b_s - math.sqrt(discriminant)) / (2 * a_s)
        
        # The valid solution for sin(t) must be in [-1, 1]
        valid_s = []
        if -1.0 <= s1 <= 1.0:
            valid_s.append(s1)
        if -1.0 <= s2 <= 1.0:
            valid_s.append(s2)
            
        if not valid_s:
            print("No valid solution for sin(t) in [-1, 1]. The chain never loses contact.")
            solution_t = None
        else:
            # We want the smallest positive time t.
            # This corresponds to the smallest positive arcsin(s).
            # Since s will be positive, we take the minimum positive s.
            s = min(s_val for s_val in valid_s if s_val >= 0)
            solution_t = math.asin(s)
            print(f"\nThe value of sin(t) when the chain loses contact is: {s:.4f}")
            print(f"The first time t when the chain loses contact is: {solution_t:.4f} seconds.")

    return solution_t

# Execute the function and store the final answer
final_time = solve_robot_chain_problem()
if final_time is not None:
    print(f'<<<{final_time:.4f}>>>')