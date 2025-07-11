import math

def calculate_critical_speed(a, b, cf, cr, m):
    """
    Calculates the critical speed for an oversteering vehicle using the linear single-track model.

    Args:
        a (float): Distance from CG to front axle (meters).
        b (float): Distance from CG to rear axle (meters).
        cf (float): Cornering stiffness of the front axle (N/rad).
        cr (float): Cornering stiffness of the rear axle (N/rad).
        m (float): Vehicle mass (kg).
    """
    print("Deriving critical speed for a single-track vehicle model.")
    print("--- Vehicle Parameters ---")
    print(f"a  (front CG distance) = {a} m")
    print(f"b  (rear CG distance)  = {b} m")
    print(f"cf (front stiffness)    = {cf} N/rad")
    print(f"cr (rear stiffness)     = {cr} N/rad")
    print(f"m  (vehicle mass)      = {m} kg\n")

    # Check for oversteering condition
    oversteer_term = a * cf - b * cr
    print(f"--- Stability Check ---")
    print(f"The vehicle is oversteering if (a * cf - b * cr) > 0.")
    print(f"({a} * {cf}) - ({b} * {cr}) = {a * cf:.2f} - {b * cr:.2f} = {oversteer_term:.2f}")

    if oversteer_term <= 0:
        print("\nThe vehicle is not oversteering. The concept of critical speed does not apply.")
        print("The vehicle is stable at all speeds according to this linear model.")
        return

    print("The vehicle is oversteering. A critical speed exists.\n")
    print("--- Critical Speed Calculation ---")
    
    wheelbase_sq = (a + b)**2
    numerator = cf * cr * wheelbase_sq
    denominator = m * oversteer_term
    
    print("The formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( (cf * cr * (a + b)^2) / (m * (a * cf - b * cr)) )\n")

    print("Substituting the given values into the equation:")
    # Print the equation with all the numbers
    equation_str = (
        f"v_crit = sqrt( ( {cf} * {cr} * ({a} + {b})^2 ) / "
        f"( {m} * ( {a} * {cf} - {b} * {cr} ) ) )"
    )
    print(equation_str)

    # Print the calculated numerator and denominator
    equation_vals = (
        f"v_crit = sqrt( {numerator:.2f} / {denominator:.2f} )"
    )
    print(equation_vals)

    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)
    v_crit_kmh = v_crit * 3.6

    print(f"v_crit = sqrt( {v_crit_sq:.2f} )")
    print(f"\nThe calculated critical speed is {v_crit:.2f} m/s, which is equivalent to {v_crit_kmh:.2f} km/h.")
    print("Above this speed, the vehicle's linear lateral dynamics become unstable.")
    
    # Final answer in the requested format
    print(f"\n<<<{v_crit:.2f}>>>")


if __name__ == '__main__':
    # --- Example Parameters ---
    # These parameters describe a slightly oversteering vehicle.
    # To be oversteering, we must satisfy: a * cf > b * cr
    
    # Distance from Center of Gravity (CG) to front axle (m)
    a_param = 1.2
    # Distance from Center of Gravity (CG) to rear axle (m)
    b_param = 1.5
    # Cornering stiffness of the front axle (N/rad)
    cf_param = 110000.0
    # Cornering stiffness of the rear axle (N/rad)
    cr_param = 85000.0
    # Vehicle mass (kg)
    m_param = 1500.0
    
    # Check condition: 1.2 * 110000 = 132000. 1.5 * 85000 = 127500. 132000 > 127500, so it is oversteering.
    
    calculate_critical_speed(a=a_param, b=b_param, cf=cf_param, cr=cr_param, m=m_param)
