import math

def calculate_critical_speed():
    """
    Calculates the critical speed for an oversteering vehicle using the linear single-track model.
    """
    # Parameters of the linear single-track model
    # a: distance from CG to front axle (m)
    # b: distance from CG to rear axle (m)
    # c_f: cornering stiffness of the front axle (N/rad)
    # c_r: cornering stiffness of the rear axle (N/rad)
    # m: vehicle mass (kg)
    # I: vehicle moment of inertia (kg*m^2) - Note: I is not needed for this specific calculation.

    # Example values for an oversteering vehicle
    # For a vehicle to be oversteering, the condition a*c_f > b*c_r must be met.
    a = 1.5   # m
    b = 1.2   # m
    c_f = 80000 # N/rad
    c_r = 90000 # N/rad
    m = 1500  # kg

    print("Vehicle Parameters:")
    print(f"a (CG to front axle) = {a} m")
    print(f"b (CG to rear axle) = {b} m")
    print(f"c_f (front cornering stiffness) = {c_f} N/rad")
    print(f"c_r (rear cornering stiffness) = {c_r} N/rad")
    print(f"m (vehicle mass) = {m} kg")
    print("-" * 30)

    # Check for the oversteering condition
    oversteer_term = a * c_f - b * c_r
    if oversteer_term <= 0:
        print("The vehicle is not oversteering (it is understeering or neutral steering).")
        print("According to the linear single-track model, it is stable at all speeds.")
        print("No critical speed exists.")
        return

    # If the vehicle is oversteering, calculate the critical speed
    print("The vehicle is oversteering. Calculating critical speed...")
    
    # Calculate numerator and denominator for the formula
    # v_crit^2 = (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r))
    numerator = c_f * c_r * (a + b)**2
    denominator = m * oversteer_term

    # Print the final equation with all the numerical values substituted
    print("\nFinal Equation with numerical values:")
    equation_str = (
        f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / "
        f"({m} * ({a} * {c_f} - {b} * {c_r})) )"
    )
    print(equation_str)

    # Calculate the critical speed
    v_crit_squared = numerator / denominator
    v_crit = math.sqrt(v_crit_squared)

    print("\nCalculation Steps:")
    print(f"Numerator = {c_f} * {c_r} * ({a+b})**2 = {numerator}")
    print(f"Denominator = {m} * ({a*c_f} - {b*c_r}) = {denominator}")
    print(f"v_crit^2 = {numerator} / {denominator} = {v_crit_squared:.2f}")
    
    print("\nResult:")
    print(f"The critical speed is {v_crit:.2f} m/s.")
    print(f"This is equivalent to {v_crit * 3.6:.2f} km/h.")
    
    # Return the final answer in the required format
    print(f"\n<<<{v_crit:.2f}>>>")


# Execute the function
calculate_critical_speed()