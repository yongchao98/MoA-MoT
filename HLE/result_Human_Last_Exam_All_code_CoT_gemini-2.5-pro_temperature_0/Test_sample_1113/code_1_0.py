import math

def solve_chain_liftoff_time():
    """
    Calculates the time when the chain first loses contact with the ground.
    """
    # Given parameters
    h = 1.0  # m, robot height
    r_cm = 25.0  # cm, arm length
    l_c = 10.0  # m, chain length
    v = 10.0  # m/s, robot speed
    beta_dot = 1.0  # rad/s, arm angular speed
    d = 20.0  # m, visible diameter of the path
    l_shadow = 10.0 * math.sqrt(3)  # m, shadow length of the path

    # --- Step 1 & 2: Calculate path and motion parameters ---
    # Radius of the circular path
    R = d / 2.0
    
    # The robot's angular speed on the path
    # omega = v / R = 10 / 10 = 1.0 rad/s
    # So, theta(t) = theta_0 + t.
    # The arm's angle beta(t) = beta_0 + beta_dot * t = pi/2 + t.
    # The problem states the robot reaches the highest point (theta=pi)
    # in the first quarter of its journey (time <= T/4 = (2*pi*R/v)/4 = pi/2).
    # t_high = pi - theta_0 <= pi/2  => theta_0 >= pi/2.
    # We assume the boundary condition to get a specific answer: theta_0 = pi/2.
    # Therefore, theta(t) = pi/2 + t.

    # Tilt angle alpha of the path
    # d * cos(alpha) = l_shadow
    cos_alpha = l_shadow / d
    alpha = math.acos(cos_alpha)
    sin_alpha = math.sin(alpha)

    # Convert arm length to meters
    r = r_cm / 100.0

    # --- Step 3 & 4: Formulate the lift-off condition ---
    # The height of the arm's tip z_tip(t) is given by:
    # z_tip(t) = R*(1-cos(theta(t)))*sin(alpha) + (h+r*cos(beta(t)))*cos(alpha) + r*sin(beta(t))*sin(theta(t))*sin(alpha)
    # With theta(t) = pi/2 + t and beta(t) = pi/2 + t, we have:
    # cos(theta) = -sin(t), sin(theta) = cos(t)
    # cos(beta) = -sin(t), sin(beta) = cos(t)
    # Substituting these gives z_tip(t) as a function of t.
    # z_tip(t) = R*(1+sin(t))*sin(alpha) + (h-r*sin(t))*cos(alpha) + r*cos^2(t)*sin(alpha)
    # The chain loses contact when z_tip(t) = l_c.
    # This can be rearranged into a quadratic equation for s = sin(t):
    # A*s^2 + B*s + C = 0
    
    A = -r * sin_alpha
    B = R * sin_alpha - r * cos_alpha
    C = (R + r) * sin_alpha + h * cos_alpha - l_c

    # --- Step 5: Solve the equation ---
    print("The problem reduces to solving for t in the equation z_tip(t) = l_c.")
    print("Let s = sin(t). The equation becomes a quadratic equation A*s^2 + B*s + C = 0.")
    print("The coefficients are:")
    print(f"A = -r*sin(alpha) = {A:.4f}")
    print(f"B = R*sin(alpha) - r*cos(alpha) = {B:.4f}")
    print(f"C = (R+r)*sin(alpha) + h*cos(alpha) - l_c = {C:.4f}")
    print("\nThe final equation to solve is:")
    print(f"{A:.4f} * s^2 + {B:.4f} * s + {C:.4f} = 0")
    
    # Calculate the discriminant
    discriminant = B**2 - 4 * A * C

    if discriminant < 0:
        print("The equation has no real solutions. The chain never loses contact.")
        return

    # Solve for s using the quadratic formula
    s1 = (-B + math.sqrt(discriminant)) / (2 * A)
    s2 = (-B - math.sqrt(discriminant)) / (2 * A)

    # --- Step 6: Calculate the time ---
    # Choose the valid solution for s (must be between -1 and 1)
    valid_s = []
    if -1.0 <= s1 <= 1.0:
        valid_s.append(s1)
    if -1.0 <= s2 <= 1.0:
        valid_s.append(s2)

    if not valid_s:
        print("No valid solution for sin(t) found.")
        return
        
    # We are looking for the first time t > 0, so we need the smallest positive s.
    s = min(s_val for s_val in valid_s if s_val > 0)
    
    # Calculate the time t
    t = math.asin(s)

    print(f"\nThe valid solution for s = sin(t) is: {s:.4f}")
    print(f"\nThe time t when the chain first loses contact with the ground is arcsin({s:.4f}).")
    print(f"t = {t:.4f} seconds")
    
    # Return the final answer in the required format
    return t

# Run the solver and print the final answer
final_time = solve_chain_liftoff_time()
if final_time is not None:
    print(f"\n<<< {final_time:.4f} >>>")
