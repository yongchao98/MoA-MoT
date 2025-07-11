import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # 1. Define the given parameters
    m = 0.1  # mass of the block in kg
    M = 10.0 # mass of the wedge in kg
    h = 2.0  # height of the wedge in meters
    theta_deg = 30.0 # angle of the wedge in degrees
    mu = 0.5 # coefficient of kinetic friction
    g = 10.0 # acceleration due to gravity in m/s^2
    # v0 is not needed as accelerations are independent of the initial velocity of the system.

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    s_theta = math.sin(theta_rad)
    c_theta = math.cos(theta_rad)

    # 2. Solve for the acceleration of the wedge (A)
    # The formula for A is derived from the equations of motion:
    # A = m*g*(mu*cos(theta)^2 - sin(theta)*cos(theta)) / (M + m*sin(theta)^2 - m*mu*sin(theta)*cos(theta))
    
    A_numerator = m * g * (mu * c_theta**2 - s_theta * c_theta)
    A_denominator = M + m * s_theta**2 - m * mu * s_theta * c_theta
    A = A_numerator / A_denominator
    
    # 3. Solve for the relative acceleration of the block (a_rel)
    # The formula for a_rel is derived from the equations of motion:
    # a_rel = g*(sin(theta) - mu*cos(theta)) - A*(mu*sin(theta) + cos(theta))
    
    term1 = g * (s_theta - mu * c_theta)
    term2 = A * (mu * s_theta + c_theta)
    a_rel = term1 - term2

    # 4. Calculate the distance the block slides (L)
    L = h / s_theta

    # 5. Calculate the time (t) using kinematics
    # L = 0.5 * a_rel * t^2  => t = sqrt(2 * L / a_rel)
    if a_rel > 0:
        time_squared = 2 * L / a_rel
        t = math.sqrt(time_squared)
    else:
        # If a_rel is zero or negative, the block does not slide down.
        t = float('inf')

    # Print the results in a clear format
    print(f"The acceleration of the wedge is A = {A:.4f} m/s^2.")
    print(f"The acceleration of the block relative to the wedge is a_rel = {a_rel:.4f} m/s^2.")
    print(f"The distance the block slides along the incline is L = {L:.4f} m.")
    print("\nThe final equation for time 't' is t = sqrt(2 * L / a_rel).")
    print("Plugging in the calculated values:")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"The final time is t = {t:.4f} s.")
    
    # Return the final answer in the specified format
    print(f"\n<<<{t}>>>")

# Execute the function to solve the problem
solve_sliding_time()