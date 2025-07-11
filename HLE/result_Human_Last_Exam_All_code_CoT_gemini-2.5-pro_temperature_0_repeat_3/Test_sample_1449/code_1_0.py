import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down a movable wedge.
    """
    # Step 1: Define given parameters
    m = 0.1  # mass of the block in kg
    M = 10.0 # mass of the wedge in kg
    theta_deg = 30.0 # angle of the wedge in degrees
    h = 2.0  # initial height of the block in m
    mu = 0.5 # coefficient of friction
    g = 10.0 # acceleration due to gravity in m/s^2

    # Convert angle to radians for use in trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 2: Calculate the relative acceleration (a_rel)
    # This formula is derived by applying Newton's second law to both the block
    # and the wedge in an inertial reference frame.
    
    # Numerator of the a_rel formula
    a_rel_numerator = g * (M + m) * (sin_theta - mu * cos_theta)
    
    # Denominator of the a_rel formula
    a_rel_denominator = M + m * sin_theta**2 - m * mu * sin_theta * cos_theta

    # The block only slides if the component of gravity along the incline
    # is greater than the maximum static friction. We assume it slides.
    # A positive numerator confirms this.
    if a_rel_numerator <= 0:
        print("The block does not slide down as the friction force is too high.")
        return

    a_rel = a_rel_numerator / a_rel_denominator

    # Step 3: Calculate the distance (L) the block slides along the incline
    L = h / sin_theta

    # Step 4: Calculate the time (t) using kinematics
    # The block starts from rest relative to the wedge, so we use L = 0.5 * a_rel * t^2
    # Solving for t gives t = sqrt(2 * L / a_rel)
    time_squared = 2 * L / a_rel
    t = math.sqrt(time_squared)

    # Output the results step-by-step as requested
    print("--- Calculation Steps ---")
    print(f"1. Distance to slide (L):")
    print(f"   L = h / sin(θ) = {h:.2f} m / sin({theta_deg}°) = {h:.2f} m / {sin_theta:.3f} = {L:.4f} m")
    
    print(f"\n2. Relative acceleration (a_rel):")
    print(f"   a_rel = Numerator / Denominator")
    print(f"   Numerator = g*(M+m)*(sin(θ) - μ*cos(θ)) = {g:.1f}*({M:.1f}+{m:.1f})*(sin({theta_deg}°) - {mu:.1f}*cos({theta_deg}°)) = {a_rel_numerator:.4f}")
    print(f"   Denominator = M + m*sin²(θ) - m*μ*sin(θ)cos(θ) = {M:.1f} + {m:.1f}*sin²({theta_deg}°) - {m:.1f}*{mu:.1f}*sin({theta_deg}°)*cos({theta_deg}°) = {a_rel_denominator:.4f}")
    print(f"   a_rel = {a_rel_numerator:.4f} / {a_rel_denominator:.4f} = {a_rel:.4f} m/s²")

    print(f"\n3. Time to slide (t):")
    print(f"   The final equation is t = sqrt(2 * L / a_rel)")
    print(f"   t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"   t = sqrt({time_squared:.4f})")
    print(f"   t = {t:.4f} s")
    
    # Final answer in the required format
    print(f"\n--- Final Answer ---")
    print(f"The exact amount of time it takes for the block to slide all the way to the bottom is {t:.4f} seconds.")
    print(f"<<<{t}>>>")

# Execute the function
solve_sliding_time()