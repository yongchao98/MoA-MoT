import math

def calculate_critical_speed(a, b, cf, cr, m):
    """
    Calculates the critical speed for a vehicle using the linear single-track model.

    Args:
        a (float): Distance from CG to front axle (m)
        b (float): Distance from CG to rear axle (m)
        cf (float): Front axle cornering stiffness (N/rad)
        cr (float): Rear axle cornering stiffness (N/rad)
        m (float): Vehicle mass (kg)
    """
    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle")
    print("="*60)
    print("The critical speed is the velocity at which the vehicle's lateral dynamics become unstable.")
    print("It is derived from the stability condition of the linear single-track model.\n")
    print("Formula: v_crit = sqrt( ((a+b)^2 * c_f * c_r) / (m * (a*c_f - b*c_r)) )\n")

    print("--- Input Parameters ---")
    print(f"a  (dist CG to front axle) = {a} m")
    print(f"b  (dist CG to rear axle)  = {b} m")
    print(f"cf (front cornering stiff) = {cf} N/rad")
    print(f"cr (rear cornering stiff)  = {cr} N/rad")
    print(f"m  (vehicle mass)          = {m} kg\n")

    # Check for oversteering condition
    front_term = a * cf
    rear_term = b * cr

    print("--- Stability Analysis ---")
    print("For instability (and a real critical speed) to exist, the vehicle must be oversteering.")
    print("Oversteer condition: a * c_f > b * c_r")
    print(f"Checking condition: {a} * {cf} > {b} * {cr}")
    print(f"Result: {front_term} > {rear_term}\n")

    if front_term > rear_term:
        print("Condition met. The vehicle is oversteering. Calculating critical speed.\n")

        # Calculations
        wheelbase_L = a + b
        numerator = (wheelbase_L**2) * cf * cr
        denominator = m * (front_term - rear_term)

        if denominator <= 0:
             print("Error: Calculation resulted in a non-positive denominator, which shouldn't happen for an oversteering vehicle.")
             return

        vcrit_squared = numerator / denominator
        vcrit_ms = math.sqrt(vcrit_squared)
        vcrit_kmh = vcrit_ms * 3.6

        # Print detailed equation steps
        print("--- Calculation Steps ---")
        print(f"v_crit = sqrt( (({a} + {b})^2 * {cf} * {cr}) / ({m} * ({a}*{cf} - {b}*{cr})) )")
        print(f"v_crit = sqrt( (({wheelbase_L})^2 * {cf*cr}) / ({m} * ({front_term} - {rear_term})) )")
        print(f"v_crit = sqrt( ({wheelbase_L**2} * {cf*cr}) / ({m} * {front_term - rear_term}) )")
        print(f"v_crit = sqrt( {numerator} / {denominator} )")
        print(f"v_crit = sqrt( {vcrit_squared:.2f} )")
        print("\n--- Result ---")
        print(f"Critical Speed (v_crit) = {vcrit_ms:.2f} m/s")
        print(f"Critical Speed (v_crit) = {vcrit_kmh:.2f} km/h")

    else:
        print("Condition not met. The vehicle is understeering or neutral steering.")
        print("According to this linear model, it is stable at all speeds.")
        print("There is no finite critical speed.")


# --- Example Parameters for an Oversteering Vehicle ---
# You can change these values to test other vehicle configurations.
a_param = 1.2   # m
b_param = 1.5   # m
cf_param = 90000.0 # N/rad
cr_param = 60000.0 # N/rad
m_param = 1500.0  # kg

# Execute the calculation
calculate_critical_speed(a_param, b_param, cf_param, cr_param, m_param)