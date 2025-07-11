import math

def derive_critical_speed():
    """
    Calculates the critical speed for an oversteering vehicle using the
    linear single-track model parameters.
    """
    # Parameters of the linear single-track model (in SI units)
    # a: distance from CG to front axle (m)
    # b: distance from CG to rear axle (m)
    # c_f: cornering stiffness of front axle (N/rad)
    # c_r: cornering stiffness of rear axle (N/rad)
    # m: vehicle mass (kg)
    # NOTE: I (moment of inertia) is not needed for the final formula.

    # Example parameters for an oversteering vehicle
    # Condition for oversteer: a*c_f > b*c_r
    a = 1.6  # m
    b = 1.2  # m
    c_f = 80000  # N/rad
    c_r = 80000  # N/rad
    m = 1500  # kg

    print("Deriving critical speed for an oversteering vehicle.")
    print("Using the following parameters (SI units):")
    print(f"a = {a} m")
    print(f"b = {b} m")
    print(f"c_f = {c_f} N/rad")
    print(f"c_r = {c_r} N/rad")
    print(f"m = {m} kg\n")

    # Check for oversteering condition
    oversteer_factor = a * c_f - b * c_r
    if oversteer_factor <= 0:
        print("The vehicle is not oversteering (a*c_f <= b*c_r).")
        print("It is stable at all speeds based on this model.")
        return

    # Calculate components for clarity
    L = a + b  # Wheelbase
    numerator = c_f * c_r * (L**2)
    denominator = m * oversteer_factor

    # Display the final equation with numerical values
    print("The formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )\n")
    print("Substituting the numerical values:")
    print(f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )")
    print(f"v_crit = sqrt( ({c_f} * {c_r} * {L**2}) / ({m} * ({a * c_f} - {b * c_r})) )")
    print(f"v_crit = sqrt( {numerator} / {denominator} )\n")

    # Calculate the critical speed
    v_crit_ms = math.sqrt(numerator / denominator)
    v_crit_kmh = v_crit_ms * 3.6

    print(f"The calculated critical speed is: {v_crit_ms:.2f} m/s")
    print(f"Which is equivalent to: {v_crit_kmh:.2f} km/h")
    
    # Final answer in the requested format
    print(f"\n<<<final_answer>>>")
    print(f"<<<{v_crit_ms:.4f}>>>")

if __name__ == "__main__":
    derive_critical_speed()