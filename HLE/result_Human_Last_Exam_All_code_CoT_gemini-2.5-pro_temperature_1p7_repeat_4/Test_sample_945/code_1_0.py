import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for a vehicle using the linear single-track model.
    The parameters are defined within the function for a sample oversteering vehicle.
    """
    # Parameters for the linear single-track model
    a = 1.2         # distance from CG to front axle (m)
    b = 1.6         # distance from CG to rear axle (m)
    cf = 120000.0   # cornering stiffness of front axle (N/rad)
    cr = 80000.0    # cornering stiffness of rear axle (N/rad)
    m = 1500.0      # vehicle mass (kg)
    
    # Note: Moment of inertia (I) and linearization speed (v) are not
    # needed for this specific critical speed calculation.

    print("Vehicle Parameters:")
    print(f"a  = {a} m")
    print(f"b  = {b} m")
    print(f"cf = {cf:.0f} N/rad")
    print(f"cr = {cr:.0f} N/rad")
    print(f"m  = {m:.0f} kg\n")

    # Step 1: Check for oversteer condition
    # A critical speed exists only for oversteering vehicles (a*cf > b*cr)
    stability_factor = a * cf - b * cr

    if stability_factor <= 0:
        print(f"The stability factor (a*cf - b*cr) is {stability_factor:.2f}.")
        print("This vehicle is not oversteering.")
        if stability_factor == 0:
            print("It is a neutral steer vehicle, which is stable at all speeds (infinite critical speed).")
        else:
            print("It is an understeer vehicle, which is stable at all speeds (no real critical speed).")
        return

    # Step 2: Calculate the critical speed using the derived formula
    # v_crit^2 = (cf * cr * (a + b)^2) / (m * (a*cf - b*cr))
    numerator = cf * cr * (a + b)**2
    denominator = m * stability_factor
    vcrit_squared = numerator / denominator
    vcrit = math.sqrt(vcrit_squared)

    # Step 3: Print the full derivation with numerical values
    print("Derivation of the Critical Speed (v_crit):")
    print("v_crit = sqrt( (cr * cf * (a + b)^2) / (m * (a*cf - b*cr)) )\n")

    print("Substituting the parameter values:")
    print(f"v_crit = sqrt( ({cr:.0f} * {cf:.0f} * ({a} + {b})^2) / ({m:.0f} * ({a}*{cf:.0f} - {b}*{cr:.0f})) )")
    print(f"v_crit = sqrt( ({cr:.0f} * {cf:.0f} * {a+b:.2f}^2) / ({m:.0f} * ({a*cf:.0f} - {b*cr:.0f})) )")
    print(f"v_crit = sqrt( {numerator:.2f} / {denominator:.2f} )")
    print(f"v_crit = sqrt( {vcrit_squared:.2f} )")
    print(f"\nResult:")
    print(f"v_crit = {vcrit:.2f} m/s")
    print(f"v_crit = {vcrit * 3.6:.2f} km/h")
    
    # This value will be used for the final answer block.
    # print(f"<<< {vcrit} >>>")

if __name__ == '__main__':
    calculate_critical_speed()
