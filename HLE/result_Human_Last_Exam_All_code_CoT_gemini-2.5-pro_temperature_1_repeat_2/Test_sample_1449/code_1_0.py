import math

def solve_sliding_time():
    """
    Calculates the time it takes for a block to slide down a movable wedge.
    """
    # --- Given Parameters ---
    m = 0.1  # Block mass in kg (100 g)
    M = 10.0 # Wedge mass in kg
    theta_deg = 30.0 # Wedge angle in degrees
    h = 2.0  # Wedge height in m
    mu = 0.5 # Coefficient of kinetic friction
    g = 10.0 # Acceleration due to gravity in m/s^2

    # --- Calculation Steps ---

    # Convert angle to radians for use in math functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 1: Calculate the distance 'L' the block slides along the incline.
    # From trigonometry, L = h / sin(theta)
    L = h / sin_theta

    # Step 2: Calculate the acceleration 'a_rel' of the block relative to the wedge.
    # The formula is derived by applying Newton's second law to both the block and the wedge.
    # a_rel = [(M+m)g*sin(theta) - (M-m)g*mu*cos(theta)] / [M + m*(sin(theta)^2 + mu*sin(theta)*cos(theta))]
    
    a_rel_numerator = (M + m) * g * sin_theta - (M - m) * g * mu * cos_theta
    a_rel_denominator = M + m * (sin_theta**2 + mu * sin_theta * cos_theta)
    a_rel = a_rel_numerator / a_rel_denominator

    # Step 3: Calculate the time 't' for the block to slide down.
    # The block starts from rest relative to the wedge.
    # Using the kinematic equation: L = 0.5 * a_rel * t^2
    # Solving for t gives: t = sqrt(2 * L / a_rel)
    time_squared = (2 * L) / a_rel
    time = math.sqrt(time_squared)

    # --- Output Results ---
    print("To find the time, we follow these steps:")
    
    print("\nStep 1: Calculate the distance 'L' the block slides.")
    print(f"L = h / sin(theta)")
    print(f"L = {h} / sin({theta_deg}°) = {h} / {sin_theta:.4f}")
    print(f"L = {L:.4f} m")

    print("\nStep 2: Calculate the relative acceleration 'a_rel'.")
    print("Numerator = (M+m)*g*sin(theta) - (M-m)*g*mu*cos(theta)")
    print(f"Numerator = ({M} + {m})*{g}*{sin_theta:.4f} - ({M} - {m})*{g}*{mu}*{cos_theta:.4f}")
    print(f"Numerator = {a_rel_numerator:.4f}")
    print("\nDenominator = M + m*(sin(theta)^2 + mu*sin(theta)*cos(theta))")
    print(f"Denominator = {M} + {m}*({sin_theta**2:.4f} + {mu}*{sin_theta:.4f}*{cos_theta:.4f})")
    print(f"Denominator = {a_rel_denominator:.4f}")
    print(f"\na_rel = Numerator / Denominator = {a_rel_numerator:.4f} / {a_rel_denominator:.4f} = {a_rel:.4f} m/s^2")

    print("\nStep 3: Calculate the final time 't'.")
    print("t = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = sqrt({time_squared:.4f})")
    print(f"Final Time, t = {time:.4f} s")

    # Return the final numerical answer for the platform
    return time

# Execute the function and capture the final answer
final_time = solve_sliding_time()
print(f"\n<<<{final_time:.4f}>>>")
