import math

def solve_slide_time():
    """
    Calculates the time for a block to slide down a movable wedge with friction.
    """
    # Given parameters
    m = 0.1  # kg (mass of the block)
    M = 10.0 # kg (mass of the wedge)
    h = 2.0  # m (initial height)
    theta_deg = 30.0 # degrees (wedge angle)
    mu = 0.5 # coefficient of kinetic friction
    g = 10.0 # m/s^2 (acceleration due to gravity)

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)

    # Pre-calculate trigonometric values
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)
    
    # 1. Derive and calculate the acceleration of the block relative to the wedge (a_rel)
    # The derivation from Newton's laws in a non-inertial frame yields:
    # a_rel = g * (sin(θ) - μ*cos(θ)) * (M + m*(sin²(θ)-cos²(θ))) / (M + m*sin²(θ) - m*μ*sin(θ)cos(θ))
    
    # Let's calculate the terms of this expression
    # Numerator part 1
    term1 = g * (sin_theta - mu * cos_theta)
    # Numerator part 2
    term2 = M + m * (sin_theta**2 - cos_theta**2)
    # Denominator
    term3 = M + m * sin_theta**2 - m * mu * sin_theta * cos_theta

    numerator = term1 * term2
    denominator = term3
    
    a_rel = numerator / denominator

    # 2. Calculate the distance 's' the block slides along the incline
    s = h / sin_theta

    # 3. Calculate the time 't' using kinematics
    # s = v_initial_rel * t + 0.5 * a_rel * t^2
    # Since v_initial_rel = 0, we have t = sqrt(2 * s / a_rel)
    time = math.sqrt(2 * s / a_rel)

    # Print the results step-by-step
    print("Step 1: Find the distance the block slides (s).")
    print(f"The incline angle is {theta_deg} degrees and the height is {h} m.")
    print(f"s = h / sin(θ) = {h} / sin({theta_deg})")
    print(f"s = {s:.4f} m\n")
    
    print("Step 2: Find the acceleration of the block relative to the wedge (a_rel).")
    print("Using Newton's laws, the formula for relative acceleration is:")
    print("a_rel = g*(sin(θ) - μ*cos(θ)) * (M + m*(sin(θ)² - cos(θ)²)) / (M + m*sin(θ)² - m*μ*sin(θ)*cos(θ))")
    print(f"With the given values:")
    print(f"m = {m} kg, M = {M} kg, θ = {theta_deg}°, μ = {mu}, g = {g} m/s²")
    print(f"a_rel = {a_rel:.4f} m/s²\n")

    print("Step 3: Calculate the time (t) it takes for the block to slide down.")
    print("Using the kinematic equation t = sqrt(2 * s / a_rel):")
    print(f"t = sqrt(2 * {s:.4f} / {a_rel:.4f})")
    print(f"t = {time:.4f} s")

solve_slide_time()
<<<3.4589>>>