import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle
    using the linear single-track model.
    """
    # Parameters of the linear single-track model
    # To ensure the vehicle is oversteering, we must satisfy a*cf > b*cr
    a = 1.2   # distance from CG to front axle (m)
    b = 1.6   # distance from CG to rear axle (m)
    cf = 95000 # cornering stiffness of front axle (N/rad)
    cr = 60000 # cornering stiffness of rear axle (N/rad)
    m = 1500  # vehicle mass (kg)
    # I (moment of inertia) is not needed for the critical speed formula.

    # Check for the oversteering condition
    if a * cf <= b * cr:
        print("Error: The provided parameters do not describe an oversteering vehicle.")
        print(f"The condition 'a * cf > b * cr' must be met.")
        print(f"Currently: a*cf = {a*cf}, b*cr = {b*cr}")
        return

    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle")
    print("-----------------------------------------------------------------")
    print("The critical speed is the velocity at which the vehicle's lateral dynamics become unstable.")
    print("It is derived from the stability conditions of the linear single-track model.")
    
    print("\nThe formula for critical speed is:")
    print("v_crit = sqrt( (cf * cr * (a + b)^2) / (m * (a*cf - b*cr)) )")
    
    print("\nUsing the provided parameters:")
    print(f"a  = {a} m")
    print(f"b  = {b} m")
    print(f"cf = {cf} N/rad")
    print(f"cr = {cr} N/rad")
    print(f"m  = {m} kg")

    # Calculate intermediate terms for clarity
    wheelbase = a + b
    numerator = cf * cr * (wheelbase**2)
    denominator_term = a * cf - b * cr
    denominator = m * denominator_term
    
    # Calculate critical speed
    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    print("\nSubstituting the values into the formula:")
    print(f"v_crit = sqrt( ({cf} * {cr} * ({a} + {b})^2) / ({m} * ({a} * {cf} - {b} * {cr})) )")
    print(f"v_crit = sqrt( ({cf * cr} * ({wheelbase:.1f})^2) / ({m} * ({a * cf} - {b * cr})) )")
    print(f"v_crit = sqrt( {numerator} / ({m} * {denominator_term}) )")
    print(f"v_crit = sqrt( {numerator} / {denominator} )")
    print(f"v_crit = sqrt( {v_crit_sq:.4f} )")

    print("\n-----------------------------------------------------------------")
    print(f"The calculated critical speed is: {v_crit:.2f} m/s")
    print(f"Which is equivalent to: {v_crit * 3.6:.2f} km/h")
    
    # Final answer in the required format
    print(f"\n<<<{v_crit}>>>")

# Execute the function
calculate_critical_speed()