import math

def solve_time_to_slide():
    """
    Calculates the time it takes for a block to slide down a moving wedge.
    """
    # --- Problem Parameters ---
    m = 0.1   # mass of the block in kg
    M = 10.0  # mass of the wedge in kg
    theta_deg = 30.0  # angle of the wedge in degrees
    h = 2.0   # height of the wedge in m
    mu = 0.5  # coefficient of friction
    g = 10.0  # acceleration due to gravity in m/s^2
    # v0 = 1.0 m/s is the initial velocity of the system and is not needed
    # for calculating the relative acceleration and time.

    # Convert angle from degrees to radians for trigonometric functions
    theta_rad = math.radians(theta_deg)
    s_t = math.sin(theta_rad)
    c_t = math.cos(theta_rad)

    print("--- Step 1: Calculate wedge acceleration A_x ---")
    # Using the formulas derived from the non-inertial frame analysis.
    # A_x = -m*g*cos(theta)*(sin(theta) - mu*cos(theta)) / (M + m*sin(theta)^2 - m*mu*sin(theta)*cos(theta))
    numerator_Ax = -m * g * c_t * (s_t - mu * c_t)
    denominator_Ax = M + m * (s_t**2) - m * mu * s_t * c_t
    Ax = numerator_Ax / denominator_Ax
    print(f"The horizontal acceleration of the wedge is A_x = {Ax:.5f} m/s^2.")

    print("\n--- Step 2: Calculate relative acceleration a_rel ---")
    # a_rel = g*(sin(theta) - mu*cos(theta)) - A_x*(cos(theta) + mu*sin(theta))
    term_gravity_friction = g * (s_t - mu * c_t)
    term_fictitious_force = Ax * (c_t + mu * s_t)
    a_rel = term_gravity_friction - term_fictitious_force
    print(f"The acceleration of the block relative to the wedge is a_rel = {a_rel:.5f} m/s^2.")

    print("\n--- Step 3: Calculate total distance L and time t ---")
    # Calculate the distance the block needs to slide
    L = h / s_t
    print(f"The distance the block slides is L = {h:.1f} / sin({theta_deg}) = {L:.5f} m.")

    # Calculate the time using the kinematic equation: L = 0.5 * a_rel * t^2
    t = math.sqrt(2 * L / a_rel)

    # Final Answer Output
    print(f"\nThe time 't' is found using the equation: t = sqrt(2 * L / a_rel)")
    print(f"t = sqrt(2 * {L:.5f} / {a_rel:.5f})")
    print(f"\nThe final time is: {t:.5f} s")

    return t

if __name__ == '__main__':
    final_time = solve_time_to_slide()
    # The final numerical answer is requested in a specific format at the end.
    # print(f"<<<{final_time:.5f}>>>") # This is for my own verification
    # I will add the final answer tag manually to the response.