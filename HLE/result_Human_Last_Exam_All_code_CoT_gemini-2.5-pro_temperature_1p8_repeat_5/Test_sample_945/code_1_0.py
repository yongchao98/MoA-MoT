import math

def calculate_critical_speed():
    """
    Calculates and prints the critical speed for an oversteering vehicle.

    The critical speed is the point at which the vehicle's linear lateral
    dynamics become unstable. This is derived from the linear single-track
    (bicycle) model of vehicle dynamics.

    Parameters for a sample oversteering vehicle:
    a  - distance from CG to front axle (m)
    b  - distance from CG to rear axle (m)
    cf - cornering stiffness of front axle (N/rad)
    cr - cornering stiffness of rear axle (N/rad)
    m  - vehicle mass (kg)
    I  - vehicle moment of inertia (kg*m^2)
    """

    # --- Vehicle Parameters ---
    # These values are chosen to satisfy the oversteer condition (a*cf > b*cr)
    a = 1.1   # m
    b = 1.6   # m
    cf = 150000 # N/rad
    cr = 80000  # N/rad
    m = 1500    # kg
    I = 3000    # kg*m^2

    # --- Oversteer Check ---
    oversteer_term = a * cf - b * cr
    if oversteer_term <= 0:
        print("The provided parameters do not describe an oversteering vehicle.")
        print(f"The term (a*cf - b*cr) must be positive, but it is {oversteer_term:.2f}.")
        return

    # --- Calculation ---
    L = a + b  # Wheelbase
    numerator = cf * cr * (L**2)
    denominator = m * oversteer_term
    
    vcrit_squared = numerator / denominator
    vcrit = math.sqrt(vcrit_squared)

    # --- Output Results ---
    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle")
    print("-" * 65)
    print("Formula:")
    print("v_crit = sqrt( (cr * cf * (a + b)^2) / (m * (a*cf - b*cr)) )")
    print("\nParameters:")
    print(f"  a  = {a} m")
    print(f"  b  = {b} m")
    print(f"  cf = {cf} N/rad")
    print(f"  cr = {cr} N/rad")
    print(f"  m  = {m} kg")
    print(f"  I  = {I} kg*m^2")
    print("\nSubstituting values into the formula:")
    
    # Print the full equation with numbers
    print(f"v_crit = sqrt( ({cr} * {cf} * ({a} + {b})^2) / ({m} * ({a}*{cf} - {b}*{cr})) )")
    
    # Print the intermediate calculation
    print(f"v_crit = sqrt( {numerator:.2f} / {denominator:.2f} )")
    print(f"v_crit = sqrt( {vcrit_squared:.2f} )")
    
    # Print the final result
    print("\nResult:")
    print(f"The critical speed is: {vcrit:.2f} m/s")

    # Final answer in the required format
    print(f"\n<<<{vcrit:.2f}>>>")

# Execute the function
calculate_critical_speed()