import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 0.1  # kg, mass of the block
    M = 10.0 # kg, mass of the wedge
    h = 2.0  # m, height of the wedge
    theta_deg = 30.0 # degrees, angle of the wedge
    mu = 0.5 # coefficient of kinetic friction
    g = 10.0 # m/s^2, acceleration due to gravity
    # v0 = 1.0 m/s, initial velocity of the system. This does not affect the
    # relative acceleration or the time to slide, as the block starts at rest
    # with respect to the wedge.

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    
    # Pre-calculate sine and cosine of the angle
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # 1. Calculate the distance the block slides, L (the length of the incline)
    L = h / sin_theta
    
    # 2. Calculate the relative acceleration of the block down the incline, a_rel.
    # This formula is derived from the equations of motion for the block and the wedge.
    # Numerator of the a_rel expression
    a_rel_num = g * (M + m) * (sin_theta - mu * cos_theta)
    
    # Denominator of the a_rel expression
    a_rel_den = M + m * sin_theta * (sin_theta - mu * cos_theta)
    
    # Check if sliding occurs (i.e., if the driving force exceeds friction)
    if a_rel_num <= 0:
        print("The block will not slide down as friction is too high or the angle is too small.")
        return

    a_rel = a_rel_num / a_rel_den
    
    # 3. Calculate the time it takes to slide down the distance L starting from rest.
    # L = 0.5 * a_rel * t^2  => t = sqrt(2 * L / a_rel)
    time_to_slide = math.sqrt(2 * L / a_rel)

    # Output the final equation with the computed numbers, as requested.
    print(f"The equation for the time (t) is: t = sqrt(2 * L / a_rel)")
    print(f"Calculated sliding distance (L): {L:.4f} m")
    print(f"Calculated relative acceleration (a_rel): {a_rel:.4f} m/s^2")
    print(f"Substituting the values into the equation:")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    
    # Final answer
    print(f"\nThe time it takes for the block to slide to the bottom is: {time_to_slide:.4f} seconds.")
    
    # Return the final numerical answer in the specified format
    return time_to_slide

# Execute the function and store the result
final_time = solve_slide_time()

# Final answer in the specified format
print(f"\n<<<{final_time}>>>")
