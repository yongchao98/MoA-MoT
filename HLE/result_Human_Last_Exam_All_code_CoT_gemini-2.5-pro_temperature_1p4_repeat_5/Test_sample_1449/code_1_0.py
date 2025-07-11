import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    
    The problem is solved by analyzing the forces in a non-inertial frame of 
    reference attached to the wedge, which yields a system of linear equations 
    for the accelerations and the normal force. Once the relative acceleration 
    of the block down the incline is found, standard kinematics is used to 
    calculate the time.
    """
    
    # --- Given values ---
    m = 0.1  # kg (100 g)
    M = 10.0 # kg
    theta_deg = 30.0 # degrees
    h = 2.0  # m
    mu = 1.0 / 2.0 # coefficient of friction
    g = 10.0 # m/s^2

    # Note: The initial velocity v0=1 m/s of the whole system is a distractor.
    # The time to slide down depends on relative accelerations, which are 
    # independent of the initial velocity of the reference frame.

    # --- Calculations ---
    
    # 1. Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # 2. Solve the system of equations for the acceleration of the wedge (A)
    # and the relative acceleration of the block (a_rel).
    #
    # From force analysis (see plan), we derive the following equations:
    # For the wedge: A = (N * (sin(theta) + mu * cos(theta))) / M
    # For the block (perpendicular to incline): N = m*g*cos(theta) - m*A*sin(theta)
    #
    # Combining these to solve for A:
    # A = (m*g*cos(theta)*(sin(theta) + mu*cos(theta))) / (M + m*sin(theta)*(sin(theta) + mu*cos(theta)))
    
    # To avoid errors, let's simplify the term m*sin(theta)*(sin(theta) + mu*cos(theta))
    # It becomes m*(sin^2(theta) + mu*sin(theta)*cos(theta))
    
    numerator_A = m * g * cos_theta * (sin_theta + mu * cos_theta)
    denominator_A = M + m * (sin_theta**2 + mu * sin_theta * cos_theta)
    A = numerator_A / denominator_A

    # 3. Now find the relative acceleration a_rel using the force equation
    # for the block parallel to the incline:
    # m*a_rel = m*g*sin(theta) - mu*N - m*A*cos(theta)
    # First, find N using the expression for A:
    N = (M * A) / (sin_theta + mu * cos_theta)
    
    # Now, solve for a_rel:
    a_rel_numerator = m * g * sin_theta - mu * N - m * A * cos_theta
    a_rel = a_rel_numerator / m
    
    # 4. Calculate the distance (L) the block slides down the incline
    L = h / sin_theta
    
    # 5. Use kinematics to find the time (t)
    # L = 0.5 * a_rel * t^2  => t = sqrt(2 * L / a_rel)
    if a_rel <= 0:
        # This would mean the block does not slide down due to friction
        time = float('inf')
    else:
        time = math.sqrt(2 * L / a_rel)

    # --- Print the final answer with the equation ---
    print("Problem parameters:")
    print(f"  Mass of block (m) = {m} kg")
    print(f"  Mass of wedge (M) = {M} kg")
    print(f"  Angle of incline (θ) = {theta_deg}°")
    print(f"  Height of wedge (h) = {h} m")
    print(f"  Coefficient of friction (μ) = {mu}")
    print(f"  Gravitational acceleration (g) = {g} m/s²\n")
    
    print("Intermediate calculations:")
    print(f"  Length of incline (L) = h / sin(θ) = {h:.1f} / {sin_theta:.3f} = {L:.4f} m")
    print(f"  Acceleration of wedge (A) = {A:.4f} m/s²")
    print(f"  Relative acceleration of block (a_rel) = {a_rel:.4f} m/s²\n")

    print("Final Calculation:")
    print(f"Time (t) = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({2 * L:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({(2 * L / a_rel):.4f})")
    print(f"t = {time:.4f} s")
    
    return time

# Execute the function to find the answer
final_time = solve_slide_time()
print(f"\n<<< {final_time} >>>")