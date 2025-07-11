import math

def solve_time():
    """
    Calculates the time for a block to slide down a moving wedge with friction.
    """
    # Given parameters
    m = 0.1  # mass of the block in kg
    M = 10.0 # mass of the wedge in kg
    theta_deg = 30.0 # angle of the wedge in degrees
    h = 2.0  # height of the wedge in m
    mu = 0.5 # coefficient of friction
    g = 10.0 # acceleration due to gravity in m/s^2

    # Convert angle to radians for math functions
    theta = math.radians(theta_deg)
    s = math.sin(theta)
    c = math.cos(theta)

    # Step 1 & 2: Solve the system of equations from the non-inertial frame analysis.
    # The derived formula for the acceleration of the wedge (A) is:
    # A = m*g*(sin(theta)*cos(theta) + mu*cos(theta)^2) / (M + m*sin(theta)^2 + mu*m*sin(theta)*cos(theta))
    numerator_A = m * g * (s * c + mu * c**2)
    denominator_A = M + m * s**2 + mu * m * s * c
    A = numerator_A / denominator_A

    # The derived formula for the acceleration of the block relative to the wedge (a) is:
    # a = g*(sin(theta) - mu*cos(theta)) - A*(cos(theta) - mu*sin(theta))
    a = g * (s - mu * c) - A * (c - mu * s)

    # Step 3: Calculate the distance L the block travels along the incline.
    # L = h / sin(theta)
    L = h / s

    # Step 4: Calculate the time t using kinematics.
    # L = 0.5 * a * t^2  => t = sqrt(2 * L / a)
    if a > 0:
        t = math.sqrt(2 * L / a)
    else:
        # This case implies the block does not slide down, which shouldn't happen
        # as tan(theta) > mu.
        print("Error: Relative acceleration is not positive, block does not slide down.")
        return

    # Output the results as requested
    print(f"The acceleration of the block relative to the wedge is a = {a:.4f} m/s^2.")
    print(f"The distance the block slides along the incline is L = {L:.4f} m.")
    print("\nThe time 't' is found using the kinematic equation L = (1/2) * a * t^2, which gives t = sqrt(2 * L / a).")
    print(f"Plugging in the values: t = sqrt(2 * {L:.4f} / {a:.4f})")
    print(f"The time it takes for the block to slide to the bottom is t = {t:.4f} s.")
    
    # Final answer format
    print(f"\n<<<{t:.4f}>>>")

solve_time()