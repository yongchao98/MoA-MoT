import math

def calculate_critical_speed():
    """
    Calculates the critical speed for an oversteering vehicle using the linear single-track model.
    """
    # Parameters of the linear single-track model
    # a: distance from CG to front axle (m)
    # b: distance from CG to rear axle (m)
    # cf: cornering stiffness of the front axle (N/rad)
    # cr: cornering stiffness of the rear axle (N/rad)
    # m: vehicle mass (kg)

    # Example values for an oversteering vehicle
    a = 1.2   # meters
    b = 1.5   # meters
    cf = 90000 # N/rad
    cr = 65000 # N/rad
    m = 1500  # kg

    # Check for the oversteering condition (a*cf > b*cr)
    # This is required for a real critical speed to exist.
    if a * cf <= b * cr:
        print("The provided parameters do not describe an oversteering vehicle.")
        print(f"Condition for oversteering (a*cf > b*cr) is not met:")
        print(f"a*cf = {a*cf:.2f}, b*cr = {b*cr:.2f}")
        return

    # Calculate the numerator of the fraction under the square root
    wheelbase = a + b
    numerator = (wheelbase**2) * cf * cr

    # Calculate the denominator of the fraction under the square root
    denominator = m * (a * cf - b * cr)

    # Calculate the critical speed
    v_crit_sq = numerator / denominator
    v_crit_ms = math.sqrt(v_crit_sq)  # result in m/s
    v_crit_kmh = v_crit_ms * 3.6       # convert to km/h

    # Print the parameters
    print("Given Vehicle Parameters:")
    print(f"a  = {a} m")
    print(f"b  = {b} m")
    print(f"cf = {cf} N/rad")
    print(f"cr = {cr} N/rad")
    print(f"m  = {m} kg")
    print("-" * 30)

    # Print the formula and the calculation steps
    print("Formula for Critical Speed (v_crit):")
    print("v_crit = sqrt( ((a+b)^2 * cf * cr) / (m * (a*cf - b*cr)) )")
    print("\nCalculation with the given values:")
    
    # Show the equation with numbers plugged in
    print(f"v_crit = sqrt( (({a} + {b})^2 * {cf} * {cr}) / ({m} * ({a}*{cf} - {b}*{cr})) )")
    
    # Show the intermediate calculation
    print(f"v_crit = sqrt( ({numerator:.2f}) / ({denominator:.2f}) )")
    
    # Show the final result
    print(f"v_crit = sqrt({v_crit_sq:.2f})")
    print(f"v_crit = {v_crit_ms:.2f} m/s  ({v_crit_kmh:.2f} km/h)")


if __name__ == "__main__":
    calculate_critical_speed()
    
<<<52.04>>>