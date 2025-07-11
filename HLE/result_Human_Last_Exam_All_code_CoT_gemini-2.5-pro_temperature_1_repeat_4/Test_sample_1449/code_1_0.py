import math

def solve_slide_time():
    """
    Calculates the time it takes for a block to slide down a movable wedge.
    """
    # Step 1: Define the given constants
    m = 0.1  # mass of the block in kg (100 g)
    M = 10.0 # mass of the wedge in kg
    h = 2.0  # height of the wedge in m
    theta_deg = 30.0 # angle of the wedge in degrees
    mu = 0.5 # coefficient of friction
    g = 10.0 # acceleration due to gravity in m/s^2

    # Step 2: Convert angle to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    sin_theta = math.sin(theta_rad)
    cos_theta = math.cos(theta_rad)

    # Step 3: Calculate the relative acceleration 'a_rel' of the block down the incline
    # This formula is derived from applying Newton's second law to the block-wedge system.
    numerator = g * (sin_theta - mu * cos_theta) * (M - m)
    denominator = M - m * sin_theta**2 + m * mu * sin_theta * cos_theta
    a_rel = numerator / denominator

    # Step 4: Calculate the distance 'L' the block slides down the incline
    L = h / sin_theta

    # Step 5: Calculate the time 't' using kinematics.
    # L = 0.5 * a_rel * t^2  =>  t = sqrt(2 * L / a_rel)
    if a_rel <= 0:
        print("The block does not slide down (or acceleration is non-positive).")
        if a_rel == 0:
            print("Force of friction and gravity component are balanced.")
        else:
            print("Friction is too high for the block to slide.")
        return

    t = math.sqrt(2 * L / a_rel)

    # Step 6: Print the results, showing the final equation with numerical values
    print("--- Intermediate Calculations ---")
    print(f"Distance to slide (L): {L:.4f} m")
    print(f"Relative acceleration (a_rel): {a_rel:.4f} m/s^2")
    print("\n--- Final Equation for Time (t) ---")
    # The problem asks to output each number in the final equation.
    print(f"t = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.4f} / {a_rel:.4f})")
    print(f"t = {t:.4f} s")
    
    return t

# Run the calculation and print the final answer
final_time = solve_slide_time()
# The final answer is returned below in the specified format.
# <<<3.4726>>>