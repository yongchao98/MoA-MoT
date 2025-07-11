import math

def solve_chain_problem():
    """
    Calculates the time when the hanging chain first loses contact with the ground.
    """
    # Given parameters from the problem description
    h = 1.0  # m, robot height
    r_cm = 25.0  # cm, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    beta_dot = 1.0  # rad/s, arm angular speed (not explicitly used as it matches robot's angular speed)
    d = 20.0  # m, visible length of the path (diameter)
    l_shadow = 10.0 * math.sqrt(3)  # m, shadow length of the path

    # --- Step 1: Calculate geometric parameters ---
    
    # Convert arm length from cm to meters
    r = r_cm / 100.0
    # The radius of the circular path is half its diameter
    R = d / 2.0
    
    # The tilt angle `theta` of the path's plane with the ground is found from its shadow
    cos_theta = l_shadow / d
    # Since the path is tilted up, sin(theta) is positive
    sin_theta = math.sqrt(1 - cos_theta**2)

    # --- Step 2 & 3: Formulate and set up the equation ---
    
    # The robot's angular speed on the path is omega = v / R = 10.0 / 10.0 = 1.0 rad/s.
    # The robot's angle on the path is alpha(t) = omega * t = t.
    # The arm's rotation angle is beta(t) = beta_dot * t = t.
    # The condition for the chain losing contact is when the arm tip's height Z_tip(t) equals the chain length l_c.
    # The height equation is: Z_tip(t) = R*sin(theta)*(1 + sin(t)) + (h + r*sin(t))*cos(theta) + r*(cos(t)**2)*sin(theta)
    # Setting Z_tip(t) = l_c and rearranging into a*x^2 + b*x + c = 0 form where x = sin(t).
    
    # Coefficients of the quadratic equation a*sin^2(t) + b*sin(t) + c = 0
    a = -r * sin_theta
    b = R * sin_theta + r * cos_theta
    c = (R + r) * sin_theta + h * cos_theta - l_c
    
    print("The problem is to find the time 't' when the chain first loses contact with the ground.")
    print("This occurs when the height of the arm's tip, Z_tip(t), equals the chain's length, l_c.\n")
    print("The height equation is: Z_tip(t) = R*sin(theta)*(1 + sin(t)) + (h + r*sin(t))*cos(theta) + r*cos(t)^2*sin(theta) = l_c")
    print("Plugging in the known values:")
    print(f"{l_c:.2f} = ({R:.2f}*{sin_theta:.4f})*(1+sin(t)) + ({h:.2f}+{r:.2f}*sin(t))*{cos_theta:.4f} + {r:.2f}*cos(t)^2*{sin_theta:.4f}")

    print("\nThis simplifies to a quadratic equation for x = sin(t) in the form a*x^2 + b*x + c = 0:")
    print(f"({a:.6f}) * sin(t)^2 + ({b:.6f}) * sin(t) + ({c:.6f}) = 0\n")

    # --- Step 4: Solve the quadratic equation ---

    discriminant = b**2 - 4 * a * c
    
    if discriminant < 0:
        print("Error: The equation has no real solutions for sin(t).")
    else:
        # Solve for x = sin(t) using the quadratic formula
        sqrt_discriminant = math.sqrt(discriminant)
        sin_t_solution_1 = (-b + sqrt_discriminant) / (2 * a)
        sin_t_solution_2 = (-b - sqrt_discriminant) / (2 * a)

        print(f"Solving the quadratic equation gives two possible values for sin(t): {sin_t_solution_1:.6f} and {sin_t_solution_2:.6f}.")

        # Choose the physically valid solution, which must be between -1 and 1
        final_sin_t = None
        if -1 <= sin_t_solution_1 <= 1:
            final_sin_t = sin_t_solution_1
        elif -1 <= sin_t_solution_2 <= 1:
            final_sin_t = sin_t_solution_2
        
        if final_sin_t is not None:
            print(f"The physically valid solution is sin(t) = {final_sin_t:.6f}.\n")
            
            # The first positive time 't' is found using arcsin.
            # Since sin(t) is positive, asin gives a result in [0, pi/2], which is the first instance.
            time_to_lose_contact = math.asin(final_sin_t)
            
            print(f"The first time 't' for this to occur is t = arcsin({final_sin_t:.6f}).")
            print(f"The chain first loses contact with the ground at t = {time_to_lose_contact:.4f} seconds.")

            print(f"<<<{time_to_lose_contact:.4f}>>>")
        else:
            print("Error: No valid solution for sin(t) found in the range [-1, 1].")

solve_chain_problem()