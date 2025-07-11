import math

def solve_chain_problem():
    """
    This function calculates the time when the chain loses contact with the ground.
    """
    # Given parameters
    h = 1.0       # m, robot height
    r = 0.25      # m, arm length
    l_c = 10.0    # m, chain length
    v = 10.0      # m/s, robot speed
    d = 20.0      # m, path diameter
    l_shadow = 10 * math.sqrt(3)  # m, shadow length

    print("--- Step 1: Analyzing the Path Geometry ---")
    # Radius of the circular path
    R = d / 2
    print(f"Path radius R = {R} m")

    # Tilt angle of the path
    cos_alpha = l_shadow / d
    alpha = math.acos(cos_alpha)
    print(f"Path tilt angle alpha = {math.degrees(alpha):.2f} degrees")
    
    # The vertical distance between the lowest and highest point of the circle is d * sin(alpha)
    # The center of the circle is at z_c = (d * sin(alpha)) / 2
    z_c = R * math.sin(alpha)
    print(f"The z-coordinate of the circle's center is {z_c:.2f} m (relative to the path's vertical midpoint).")


    print("\n--- Step 2: Describing the Robot's Motion ---")
    # Angular velocity of the robot on the path
    omega_robot = v / R
    print(f"Robot's angular velocity on the path omega = {omega_robot:.2f} rad/s")
    # With omega=1, the angle is equal to time t. We assume the robot starts at an angle
    # that leads it to the highest point at t=pi/2, which simplifies the motion equations.
    # The z-coordinate of the robot's base can be expressed as:
    # z_base(t) = z_c + z_c * sin(t) = 5 * (1 + sin(t))
    print("The height of the robot's base is z_base(t) = 5 * (1 + sin(t))")


    print("\n--- Step 3: Finding the Height of the Arm's Tip ---")
    # The z-component of the normal vector to the plane (pointing "up") is sqrt(3)/2
    n_z = math.sqrt(3) / 2
    # The height of the robot's top is z_top(t) = z_base(t) + h * n_z
    print("The height of the robot's top is z_top(t) = z_base(t) + h * n_z = 5*(1+sin(t)) + sqrt(3)/2")

    # The z-component of the arm's vector A(t) is z_A(t) = r * (1/2 * cos^2(t) + sqrt(3)/2 * sin(t))
    print("The height contribution from the arm is z_arm(t) = 0.25 * (0.5*cos^2(t) + (sqrt(3)/2)*sin(t))")
    # So, the total height of the arm's tip is z_tip(t) = z_top(t) + z_A(t)
    

    print("\n--- Step 4: Solving for the Time of Detachment ---")
    print(f"The chain loses contact when the tip's height z_tip(t) equals the chain length l_c = {l_c} m.")
    print("The equation is: 5*(1+sin(t)) + sqrt(3)/2 + 0.25*(0.5*(1-sin^2(t)) + (sqrt(3)/2)*sin(t)) = 10")
    print("Let s = sin(t). After simplification, this yields a quadratic equation: a*s^2 + b*s + c = 0")
    
    # Coefficients of the quadratic equation s^2 - (40+sqrt(3))s + (39-4sqrt(3)) = 0
    a = 1.0
    b = -(40 + math.sqrt(3))
    c = 39 - 4 * math.sqrt(3)

    print("\nThe equation is of the form a*s^2 + b*s + c = 0, where:")
    print(f"a = {a}")
    print(f"b = -(40 + \u221A3) = {b:.4f}")
    print(f"c = 39 - 4*\u221A3 = {c:.4f}")

    print("\nTherefore, the final equation in terms of s = sin(t) is:")
    print(f"{a} * s\u00B2 - {40 + math.sqrt(3):.4f} * s + {39 - 4 * math.sqrt(3):.4f} = 0")
    
    # Solve the quadratic equation for s
    discriminant = b**2 - 4 * a * c
    s1 = (-b + math.sqrt(discriminant)) / (2 * a)
    s2 = (-b - math.sqrt(discriminant)) / (2 * a)
    
    print(f"\nThe two potential solutions for s are: {s1:.4f} and {s2:.4f}")

    # The valid solution for s must be between -1 and 1
    if -1 <= s2 <= 1:
        s_valid = s2
    else:
        s_valid = s1

    print(f"Since s = sin(t), its value must be between -1 and 1. The valid solution is s = {s_valid:.4f}")
    
    # Calculate time t from s
    t = math.asin(s_valid)

    print("\n--- Final Result ---")
    print(f"The time t when the chain first loses contact is arcsin({s_valid:.4f})")
    print(f"t = {t:.3f} seconds.")
    
    return t

# Run the solver
final_time = solve_chain_problem()
# Return the final answer in the specified format
# The value is ~0.8998, which we can round.
print(f'<<<{final_time:.3f}>>>')
