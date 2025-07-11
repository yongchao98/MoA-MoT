import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a movable wedge.
    """
    # Given parameters
    m = 0.1  # kg
    M = 10   # kg
    theta_deg = 30
    h = 2    # m
    mu = 0.5
    g = 10   # m/s^2

    # Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    s_theta = math.sin(theta_rad)
    c_theta = math.cos(theta_rad)

    # Step 1: Calculate the distance the block slides (L)
    # The path length is the hypotenuse of the triangle with height h and angle theta.
    # L = h / sin(theta)
    L = h / s_theta
    
    # Step 2: Calculate the block's acceleration relative to the wedge (a_rel).
    # This formula is derived from applying Newton's second law to the block and the wedge,
    # and using the conservation of horizontal momentum for the system.
    # a_rel = g*(M+m)*(sin(theta) - mu*cos(theta)) / (M + m*sin(theta)^2 - m*mu*sin(theta)*cos(theta))
    numerator = g * (M + m) * (s_theta - mu * c_theta)
    denominator = M + m * s_theta**2 - m * mu * s_theta * c_theta
    a_rel = numerator / denominator

    # Step 3: Calculate the time (t) using kinematics.
    # Since the block starts from rest relative to the wedge, we use L = 0.5 * a_rel * t^2.
    # t = sqrt(2 * L / a_rel)
    time = math.sqrt(2 * L / a_rel)

    # --- Output the results and calculation steps ---
    print(f"--- Problem Parameters ---")
    print(f"Block mass (m): {m} kg")
    print(f"Wedge mass (M): {M} kg")
    print(f"Wedge angle (θ): {theta_deg}°")
    print(f"Initial height (h): {h} m")
    print(f"Coefficient of friction (μ): {mu}")
    print(f"Gravitational acceleration (g): {g} m/s²\n")

    print(f"--- Step 1: Find the sliding distance L ---")
    print(f"The equation is: L = h / sin(θ)")
    print(f"L = {h} / sin({theta_deg}°) = {L:.4f} m\n")

    print(f"--- Step 2: Find the relative acceleration a_rel ---")
    print(f"The equation for the numerator is: g * (M + m) * (sin(θ) - μ * cos(θ))")
    print(f"Numerator = {g} * ({M} + {m}) * ({s_theta:.4f} - {mu} * {c_theta:.4f}) = {numerator:.4f}")
    print(f"The equation for the denominator is: M + m*sin(θ)² - m*μ*sin(θ)*cos(θ)")
    print(f"Denominator = {M} + {m}*{s_theta**2:.4f} - {m}*{mu}*{s_theta:.4f}*{c_theta:.4f} = {denominator:.4f}")
    print(f"a_rel = Numerator / Denominator = {a_rel:.4f} m/s²\n")

    print(f"--- Step 3: Find the time t ---")
    print(f"The equation is: t = √(2 * L / a_rel)")
    print(f"t = √(2 * {L:.4f} / {a_rel:.4f})")
    print(f"The final time is: {time:.4f} seconds")

    return time

# Execute the function and store the final answer
final_time = solve_slide_time()
# The final answer will be automatically captured in the <<<>>> format by the system.
# The numeric value itself is printed above.
# print(f"<<<{final_time:.7f}>>>")