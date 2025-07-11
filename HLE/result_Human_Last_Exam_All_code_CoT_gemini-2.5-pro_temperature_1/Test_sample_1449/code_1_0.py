import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Step 1: Define the given physical constants.
    m = 100 / 1000  # mass of the block in kg
    M = 10.0        # mass of the wedge in kg
    theta_deg = 30.0  # angle of the wedge in degrees
    h = 2.0         # initial height of the block in meters
    mu = 1.0 / 2.0    # coefficient of kinetic friction
    g = 10.0        # acceleration due to gravity in m/s^2

    # Step 2: Convert the angle from degrees to radians for use in trigonometric functions.
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 3: Calculate the relative acceleration 'a_rel' of the block down the incline.
    # This formula is derived from Newton's laws in a non-inertial reference frame
    # and the conservation of horizontal momentum.
    # a_rel = g * (sin(θ) - μ*cos(θ)) / (1 - (m/(m+M)) * cos(θ) * (μ*sin(θ) + cos(θ)))

    # Calculate the numerator of the expression for a_rel
    numerator_a_rel = g * (sin_theta - mu * cos_theta)

    # Calculate the denominator of the expression for a_rel
    denominator_a_rel = 1 - (m / (m + M)) * cos_theta * (mu * sin_theta + cos_theta)

    # Calculate the final relative acceleration
    a_rel = numerator_a_rel / denominator_a_rel

    # Step 4: Calculate the distance 'L' the block slides along the incline.
    L = h / sin_theta

    # Step 5: Calculate the time 't' to slide down distance 'L' starting from relative rest.
    # The kinematic equation is L = 0.5 * a_rel * t^2.
    # We solve for t: t = sqrt(2 * L / a_rel).
    t_squared = (2 * L) / a_rel
    if t_squared < 0:
        print("Error: Calculation results in a negative value under the square root.")
        print("This implies the block does not slide down (friction is too high).")
        return
        
    t = math.sqrt(t_squared)

    # Step 6: Print the final equation with the calculated numbers as requested.
    print("The final equation for time 't' is: t = sqrt(2 * L / a_rel)")
    print("---------------------------------------------------------")
    print(f"Calculated distance to slide, L = h / sin(θ) = {h:.1f} / {sin_theta:.3f} = {L:.4f} m")
    print(f"Calculated relative acceleration, a_rel = {a_rel:.4f} m/s^2")
    print("---------------------------------------------------------")
    print(f"Plugging in the numbers: t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"The time it takes for the block to slide down is: {t:.4f} s")
    
    # Final answer in the required format
    print(f"\n<<<{t}>>>")

solve_slide_time()