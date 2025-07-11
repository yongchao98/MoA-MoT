import math

def calculate_critical_speed():
    """
    Calculates the critical speed for an oversteering vehicle based on the
    linear single-track model.
    """
    # Parameters of the linear single-track model
    # a: distance from CG to front axle (m)
    # b: distance from CG to rear axle (m)
    # c_f: cornering stiffness of the front axle (N/rad)
    # c_r: cornering stiffness of the rear axle (N/rad)
    # m: vehicle mass (kg)
    # I: vehicle moment of inertia (kg*m^2) - Note: not needed for critical speed
    
    # --- Example values for an oversteering vehicle ---
    a = 1.1   # meters
    b = 1.6   # meters
    c_f = 90000 # N/rad
    c_r = 50000 # N/rad
    m = 1500  # kg
    # I = 2500, not used in the final formula but is part of the model.

    print("--- Vehicle Parameters ---")
    print(f"a (front CG distance): {a} m")
    print(f"b (rear CG distance): {b} m")
    print(f"c_f (front cornering stiffness): {c_f} N/rad")
    print(f"c_r (rear cornering stiffness): {c_r} N/rad")
    print(f"m (vehicle mass): {m} kg\n")

    # Check for oversteering condition: a*c_f > b*c_r
    # If this is not met, the vehicle is neutral or understeering,
    # and it does not have a critical speed (it's always stable).
    front_factor = a * c_f
    rear_factor = b * c_r

    print("--- Stability Analysis ---")
    print("An oversteering vehicle becomes unstable at a certain 'critical speed'.")
    print("A vehicle is oversteering if a*c_f > b*c_r.")
    print(f"Checking condition: {a} * {c_f} > {b} * {c_r}")
    print(f"Result: {front_factor} > {rear_factor}\n")

    if front_factor <= rear_factor:
        print("Result: The vehicle is not oversteering. It is stable at all speeds according to the linear model.")
        return

    # Calculate critical speed for the oversteering vehicle
    # v_crit^2 = (c_f * c_r * (a+b)^2) / (m * (a*c_f - b*c_r))
    wheelbase = a + b
    numerator = c_f * c_r * (wheelbase ** 2)
    denominator_term = front_factor - rear_factor
    denominator = m * denominator_term

    print("--- Critical Speed Calculation ---")
    print("The formula for the square of the critical speed is:")
    print("v_crit^2 = (c_f * c_r * (a+b)^2) / (m * (a*c_f - b*c_r))\n")
    
    print("Plugging in the numbers:")
    # Using f-strings to display the equation with values
    print(f"v_crit^2 = ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a}*{c_f} - {b}*{c_r}))")
    print(f"v_crit^2 = ({c_f} * {c_r} * {wheelbase:.2f}^2) / ({m} * ({front_factor} - {rear_factor}))")
    print(f"v_crit^2 = ({numerator:.4e}) / ({m} * {denominator_term})")
    print(f"v_crit^2 = {numerator:.4e} / {denominator:.4e}\n")

    v_crit_squared = numerator / denominator
    v_crit_ms = math.sqrt(v_crit_squared)
    v_crit_kmh = v_crit_ms * 3.6
    
    print(f"v_crit^2 = {v_crit_squared:.2f} m^2/s^2")
    print(f"v_crit = sqrt({v_crit_squared:.2f})")
    
    print("\n--- Final Result ---")
    print(f"Critical Speed (v_crit): {v_crit_ms:.2f} m/s")
    print(f"Critical Speed (v_crit): {v_crit_kmh:.2f} km/h")
    
    # Final answer for the script
    # This will be the value in m/s, which is the standard unit.
    global final_answer
    final_answer = v_crit_ms


# Run the calculation and store the result
final_answer = None
calculate_critical_speed()
# The required output format.
if final_answer is not None:
    print(f'<<<{final_answer:.2f}>>>')
