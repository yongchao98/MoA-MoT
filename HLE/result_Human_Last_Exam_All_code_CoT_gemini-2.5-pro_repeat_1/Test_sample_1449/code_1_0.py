import math

def solve_time():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # Given parameters
    m = 100 / 1000  # Convert grams to kg
    M = 10.0        # kg
    g = 10.0        # m/s^2
    h = 2.0         # m
    theta_deg = 30.0
    mu = 0.5

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # --- Step 1: Calculate the acceleration of the block relative to the wedge (a_rel) ---
    # This formula is derived by conserving horizontal momentum and applying Newton's laws.
    # a_rel = (g * (M + m) * (sin(theta) - mu * cos(theta))) / 
    #         (M + m * (sin(theta)**2 - mu * sin(theta) * cos(theta)))
    
    a_rel_numerator = (M + m) * g * (sin_theta - mu * cos_theta)
    a_rel_denominator = M + m * (sin_theta**2 - mu * sin_theta * cos_theta)
    
    if a_rel_denominator == 0:
        print("Error: Division by zero in acceleration calculation.")
        return

    a_rel = a_rel_numerator / a_rel_denominator

    if a_rel <= 0:
        print("The block does not slide down as the friction is too high or the slope is too shallow.")
        return

    # --- Step 2: Calculate the distance the block slides along the wedge (L) ---
    L = h / sin_theta

    # --- Step 3: Calculate the time taken (t) using kinematics ---
    # The block starts from rest relative to the wedge, so L = 0.5 * a_rel * t^2
    # Solving for t gives: t = sqrt(2 * L / a_rel)
    time_squared = 2 * L / a_rel
    t = math.sqrt(time_squared)

    # --- Output the results as requested ---
    print("The final time is calculated using the equation: t = sqrt(2 * L / a_rel)")
    print(f"1. The distance the block slides, L = h / sin({theta_deg}) = {h:.1f} / {sin_theta:.1f} = {L:.4f} m.")
    print(f"2. The relative acceleration of the block, a_rel = {a_rel:.4f} m/s^2.")
    print(f"Plugging the values into the equation:")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = {t:.4f} s")
    
    # Return the final answer in the required format
    print(f"\n<<< {t:.4f} >>>")

solve_time()