import math

def calculate_critical_speed(a, b, cf, cr, m, I):
    """
    Calculates the critical speed for an oversteering vehicle.

    Parameters:
    a (float): distance from CG to front axle (m)
    b (float): distance from CG to rear axle (m)
    cf (float): cornering stiffness of the front axle (N/rad)
    cr (float): cornering stiffness of the rear axle (N/rad)
    m (float): vehicle mass (kg)
    I (float): vehicle moment of inertia (kg*m^2) - Note: Not needed for this formula but included for completeness.
    """

    # Condition for oversteer: b*cr > a*cf
    oversteer_check = b * cr - a * cf
    if oversteer_check <= 0:
        print("The vehicle is not oversteering based on the provided parameters.")
        print("The condition b*cr > a*cf is not met.")
        print(f"{b} * {cr} (which is {b*cr}) is not greater than {a} * {cf} (which is {a*cf}).")
        return

    # Numerator of the fraction under the square root
    numerator = (a + b)**2 * cf * cr
    # Denominator of the fraction under the square root
    denominator = m * oversteer_check

    # Calculate critical speed
    vcrit_squared = numerator / denominator
    vcrit = math.sqrt(vcrit_squared)

    # Print the equation with the values
    print("Derivation for critical speed (vcrit) of an oversteering vehicle:")
    print("vcrit = sqrt( ((a + b)^2 * cf * cr) / (m * (b * cr - a * cf)) )")
    print(f"vcrit = sqrt( (({a} + {b})^2 * {cf} * {cr}) / ({m} * ({b} * {cr} - {a} * {cf})) )")
    print(f"vcrit = sqrt( ({numerator}) / ({denominator}) )")
    print(f"vcrit = sqrt({vcrit_squared})")
    print(f"Final Critical Speed: {vcrit} m/s")
    
    return vcrit

# --- Example Parameters ---
# Feel free to change these values to test different scenarios
a_param = 1.1  # meters
b_param = 1.4  # meters
cf_param = 80000  # N/rad
cr_param = 120000 # N/rad
m_param = 1500  # kg
I_param = 3000   # kg*m^2 (not used in the final formula)

# --- Calculation and Output ---
# This ensures that b*cr (168000) > a*cf (88000), so the vehicle is oversteering
final_vcrit = calculate_critical_speed(a_param, b_param, cf_param, cr_param, m_param, I_param)

# The following line is to fulfill the final output format requirement.
# For example use: <<<27.95>>>
print(f'<<<{final_vcrit}>>>')