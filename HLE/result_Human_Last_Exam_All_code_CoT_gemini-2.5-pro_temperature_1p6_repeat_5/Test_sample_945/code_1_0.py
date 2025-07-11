import math

def calculate_critical_speed(a, b, c_f, c_r, m):
    """
    Calculates the critical speed for a vehicle using the linear single-track model.
    
    Parameters:
    a   (float): distance from CG to front axle (m)
    b   (float): distance from CG to rear axle (m)
    c_f (float): cornering stiffness of the front axle (N/rad)
    c_r (float): cornering stiffness of the rear axle (N/rad)
    m   (float): vehicle mass (kg)
    """

    print("--- Vehicle Parameters ---")
    print(f"a (CG to front axle): {a} m")
    print(f"b (CG to rear axle):  {b} m")
    print(f"c_f (front stiffness):  {c_f} N/rad")
    print(f"c_r (rear stiffness):   {c_r} N/rad")
    print(f"m (mass):             {m} kg\n")
    
    # Check for oversteer condition
    front_term = a * c_f
    rear_term = b * c_r
    
    print("--- Stability Condition Check ---")
    print(f"Condition for oversteer: a * c_f > b * c_r")
    print(f"Checking values: {front_term:.2f} > {rear_term:.2f}")

    if front_term <= rear_term:
        print("\nResult: The vehicle is not oversteering. It is understeering or neutral steering.")
        print("Therefore, it does not have a critical speed and is stable at all speeds in this linear model.")
        return

    print("Result: The vehicle is oversteering and has a critical speed.\n")
    
    # --- Calculation ---
    print("--- Critical Speed Calculation ---")
    print("Formula: v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )\n")

    # Substitute values into the formula
    numerator = c_f * c_r * (a + b)**2
    oversteer_factor = front_term - rear_term
    denominator = m * oversteer_factor
    
    print("Step 1: Substitute parameters into the formula")
    print(f"v_crit = sqrt( ({c_r} * {c_f} * ({a} + {b})^2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )")
    
    print("\nStep 2: Calculate terms in the equation")
    print(f"v_crit = sqrt( ({c_r} * {c_f} * {a+b:.2f}^2) / ({m} * ({front_term:.2f} - {rear_term:.2f})) )")
    print(f"v_crit = sqrt( ({numerator:.2f}) / ({m} * {oversteer_factor:.2f}) )")
    print(f"v_crit = sqrt( {numerator:.2f} / {denominator:.2f} )")

    if denominator <= 0:
        print("\nError: Calculation resulted in a non-positive denominator. Cannot compute a real critical speed.")
        return

    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    print(f"v_crit = sqrt({v_crit_sq:.4f})")
    
    print("\n--- Final Result ---")
    print(f"Critical Speed (v_crit): {v_crit:.4f} m/s")
    print(f"Equivalent to: {v_crit * 3.6:.4f} km/h")
    
    # Final answer in the specified format
    print(f"\n<<<{v_crit:.4f}>>>")

# --- Define Vehicle Parameters Here ---
# You can change these values to analyze a different vehicle.
# These parameters represent an oversteering vehicle.
a_param = 1.1     # m
b_param = 1.6     # m
c_f_param = 95000 # N/rad
c_r_param = 60000 # N/rad
m_param = 1500    # kg

# Run the calculation
calculate_critical_speed(a_param, b_param, c_f_param, c_r_param, m_param)