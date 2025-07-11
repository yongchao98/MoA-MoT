import math

def solve_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    lc = 10.0  # m, chain length
    d = 20.0  # m, visible path length (diameter)
    l_shadow = 10.0 * math.sqrt(3)  # m, shadow length
    v = 10.0  # m/s, robot speed
    beta_dot = 1.0  # rad/s, arm angular speed

    # Step 1: Determine the geometry of the path
    R = d / 2.0  # Radius of the circular path
    # d * cos(theta) = l_shadow
    cos_theta = l_shadow / d
    # Clamp the value to handle potential floating point inaccuracies
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.acos(cos_theta)  # Tilt angle in radians
    sin_theta = math.sin(theta)

    # Step 2: Define kinematics
    # Robot angular speed along the path: omega = v / R
    omega = v / R
    # Since omega = 10/10 = 1 rad/s, the path angle alpha(t) = t.
    # Since beta_dot = 1 rad/s, the arm angle beta(t) = t.
    
    # Step 3: Formulate the height equation
    # The equation for the height of the arm's tip (Z_tip) is:
    # Z_tip(t) = R*sin(theta)*(1+sin(t)) + h*cos(theta) + r*(cos^2(t)*sin(theta) + sin(t)*cos(theta))
    # The chain loses contact when Z_tip(t) = lc.
    # We get the equation:
    # R*sin(theta)*(1+sin(t)) + h*cos(theta) + r*(cos^2(t)*sin(theta) + sin(t)*cos(theta)) = lc
    #
    # Substitute cos^2(t) = 1 - sin^2(t) and let x = sin(t).
    # This leads to a quadratic equation: a*x^2 + b*x + c = 0
    
    # Coefficients for the quadratic equation a*x^2 + b*x + c = 0 where x = sin(t)
    a = -r * sin_theta
    b = R * sin_theta + r * cos_theta
    c = R * sin_theta + h * cos_theta + r * sin_theta - lc

    print("The problem is solved by finding the time 't' that satisfies the height equation:")
    print("Height_of_arm_tip(t) = Length_of_chain")
    print("\nThe equation with the given values is:")
    print(f"({R:.1f}*sin({math.degrees(theta):.0f}°))*(1 + sin(t)) + ({h:.1f}*cos({math.degrees(theta):.0f}°)) + "
          f"{r:.2f}*(cos²(t)*sin({math.degrees(theta):.0f}°) + sin(t)*cos({math.degrees(theta):.0f}°)) = {lc:.1f}")
    
    # Substitute numerical values
    term1_val = R * sin_theta
    term2_val = h * cos_theta
    term3_r_sin = r * sin_theta
    term3_r_cos = r * cos_theta

    print("\nSubstituting the numerical values for the trigonometric functions:")
    print(f"{term1_val:.4f} * (1 + sin(t)) + {term2_val:.4f} + {r:.2f} * (cos²(t)*{sin_theta:.4f} + sin(t)*{cos_theta:.4f}) = {lc:.1f}")
    
    print("\nThis simplifies to a quadratic equation for x = sin(t) of the form ax² + bx + c = 0:")
    print(f"({a:.4f})x² + ({b:.4f})x + ({c:.4f}) = 0")

    # Step 4: Solve the quadratic equation
    discriminant = b**2 - 4*a*c

    if discriminant < 0:
        print("\nNo real solution for time exists.")
        return

    # Calculate the two possible solutions for x = sin(t)
    sqrt_discriminant = math.sqrt(discriminant)
    x1 = (-b + sqrt_discriminant) / (2*a)
    x2 = (-b - sqrt_discriminant) / (2*a)

    # Find the valid solution for sin(t) which must be between -1 and 1
    valid_sin_t = None
    if -1.0 <= x1 <= 1.0:
        valid_sin_t = x1
    elif -1.0 <= x2 <= 1.0:
        valid_sin_t = x2

    if valid_sin_t is None:
        print("\nNo valid solution for sin(t) in the range [-1, 1].")
        return
        
    # Calculate the first time 't'
    t = math.asin(valid_sin_t)

    # We need the first positive time. If t is negative, the next positive solution is pi - t.
    # However, given the setup, the first solution should be positive.
    if t < 0:
        # This case is unlikely for this problem's physical setup
        t = math.pi - t

    print(f"\nThe valid solution for sin(t) is: {valid_sin_t:.4f}")
    print(f"\nThe first time 't' when the chain loses contact with the ground is arcsin({valid_sin_t:.4f}).")
    print(f"\nt = {t:.4f} seconds")
    
    # Return the final answer in the required format
    return t

# Run the solver and capture the final numerical answer
final_time = solve_chain_problem()
print(f"\n<<<{final_time:.4f}>>>")
