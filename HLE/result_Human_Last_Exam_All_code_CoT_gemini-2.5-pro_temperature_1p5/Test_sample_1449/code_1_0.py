import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down an accelerating wedge with friction.
    """
    # 1. Define constants from the problem statement
    m = 0.1  # mass of the block in kg (100 g)
    M = 10.0 # mass of the wedge in kg
    h = 2.0  # height of the wedge in meters
    theta_deg = 30.0 # angle of the wedge in degrees
    mu = 0.5 # coefficient of friction
    g = 10.0 # acceleration due to gravity in m/s^2
    # v0 = 1.0 m/s is the initial velocity of the entire system and does not affect the relative motion.

    # 2. Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    
    # 3. Calculate the sliding distance 'd'
    d = h / math.sin(theta_rad)
    print(f"The block will slide a distance d = h / sin({theta_deg}°) = {h:.1f} / {math.sin(theta_rad):.3f} = {d:.4f} m.")

    # 4. Calculate the relative acceleration 'a_rel' using the derived formula
    s = math.sin(theta_rad)
    c = math.cos(theta_rad)
    c2 = math.cos(2 * theta_rad)
    
    # Numerator of the acceleration formula
    num_term1 = M * (s - mu * c)
    num_term2 = m * c2 * (s + mu * c)
    numerator = g * (num_term1 - num_term2)
    
    # Denominator of the acceleration formula
    den_term1 = m * (s**2 + mu * s * c)
    denominator = M + den_term1
    
    a_rel = numerator / denominator
    print(f"The relative acceleration of the block down the wedge, a_rel, is {a_rel:.4f} m/s^2.")

    # 5. Calculate the time 't' using kinematics
    # d = 0.5 * a_rel * t^2  => t = sqrt(2 * d / a_rel)
    t = math.sqrt(2 * d / a_rel)
    
    print("\nCalculating the final time 't':")
    print(f"t = sqrt(2 * d / a_rel)")
    print(f"t = sqrt(2 * {d:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({2 * d / a_rel:.4f})")
    print(f"The exact time it takes for the block to slide down is {t:.4f} seconds.")

solve_sliding_time()
<<<3.5912>>>