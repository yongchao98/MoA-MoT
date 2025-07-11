import math

def solve_robot_chain_problem():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Parameters from the problem statement
    h = 1.0  # m, robot height
    r = 0.25  # m, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    d = 20.0  # m, visible diameter of the path
    l_shadow = 10 * math.sqrt(3)  # m, shadow length

    # Step 1: Determine geometric parameters
    R = d / 2.0  # Radius of the circular path
    # The shadow length l_shadow is the minor axis of the elliptical projection.
    # l_shadow = d * cos(alpha) = 2 * R * cos(alpha)
    cos_alpha = l_shadow / d
    sin_alpha = math.sqrt(1 - cos_alpha**2)
    
    # Step 2 & 3: Formulate the height equation and analyze motion
    # The height of the arm tip z_tip(t) is given by:
    # z_tip = R*sin(a)*(1+sin(phi)) + h*cos(a) + r*(cos(t)*cos(phi)*sin(a) + sin(t)*cos(a))
    # Motion: phi(t) = phi(0) + (v/R)*t = phi(0) + t. arm_angle(t) = t.
    # The chain loses contact when z_tip(t) = l_c.
    
    # Step 4: Resolve the starting position phi(0)
    # The problem is only solvable for a specific phi(0). We assume the problem implies
    # the boundary case of the "first quarter" condition, which means phi(0) = 0.
    # With phi(0)=0, phi(t)=t. The equation becomes:
    # R*sin(a)*(1+sin(t)) + h*cos(a) + r*(cos(t)*cos(t)*sin(a) + sin(t)*cos(a)) = l_c

    print("The final equation to solve for t is:")
    print(f"{R}*sin(alpha)*(1+sin(t)) + {h}*cos(alpha) + {r}*(cos(t)^2*sin(alpha) + sin(t)*cos(alpha)) = {l_c}")
    print("Substituting the known values (alpha = arccos(sqrt(3)/2) = pi/6):")
    print(f"{R}*{sin_alpha:.1f}*(1+sin(t)) + {h}*{cos_alpha:.4f} + {r}*(cos(t)^2*{sin_alpha:.1f} + sin(t)*{cos_alpha:.4f}) = {l_c}")
    print("\nThis simplifies to a quadratic equation for x = sin(t).")

    # Step 5: Solve for t
    # Let x = sin(t). Then cos^2(t) = 1 - x^2.
    # The equation is:
    # R*sin_alpha*(1+x) + h*cos_alpha + r*((1-x^2)*sin_alpha + x*cos_alpha) = l_c
    # Rearranging into ax^2 + bx + c = 0 form:
    # x^2*(-r*sin_alpha) + x*(R*sin_alpha + r*cos_alpha) + (R*sin_alpha + h*cos_alpha + r*sin_alpha - l_c) = 0
    
    a = -r * sin_alpha
    b = R * sin_alpha + r * cos_alpha
    c = R * sin_alpha + h * cos_alpha + r * sin_alpha - l_c

    # Solve the quadratic equation for x = sin(t)
    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("No real solution for time exists.")
        return

    # We need the smaller positive solution for sin(t) for the first time t>0.
    # The quadratic formula is (-b +/- sqrt(discriminant)) / (2a)
    sin_t_1 = (-b + math.sqrt(discriminant)) / (2 * a)
    sin_t_2 = (-b - math.sqrt(discriminant)) / (2 * a)

    # Choose the physically meaningful solution for sin(t) which must be between -1 and 1
    valid_sin_t = None
    if -1 <= sin_t_1 <= 1:
        valid_sin_t = sin_t_1
    if -1 <= sin_t_2 <= 1:
        # If both are valid, we usually take the one that gives the smallest positive t.
        # In this case, sin_t_2 is the larger value.
        valid_sin_t = sin_t_2
        
    if valid_sin_t is None:
        print("No valid solution for sin(t) in the range [-1, 1].")
        return

    # Calculate time t
    time = math.asin(valid_sin_t)

    print("\nSolving the quadratic equation gives sin(t) =", round(valid_sin_t, 4))
    print(f"The time t is arcsin({round(valid_sin_t, 4)}) seconds.")
    print("\nFinal Answer:")
    print(time)

solve_robot_chain_problem()