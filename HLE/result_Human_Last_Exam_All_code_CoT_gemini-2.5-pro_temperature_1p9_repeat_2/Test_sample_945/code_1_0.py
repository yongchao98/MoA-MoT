import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle based on the
    linear single-track model.
    """
    # Parameters of the linear single-track model
    # Using example values for a sporty passenger car configured for oversteer
    a = 1.2   # distance between CG and front axle (m)
    b = 1.5   # distance between CG and rear axle (m)
    cf = 90000 # cornering stiffness of the front axle (N/rad)
    cr = 65000 # cornering stiffness of the rear axle (N/rad)
    m = 1600  # vehicle mass (kg)
    # Note: The vehicle moment of inertia 'I' is not needed for the final formula
    # as it cancels out during the derivation of the static stability limit.

    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("---------------------------------------------------------")
    print("The critical speed (v_crit) is the speed at which the linear lateral dynamics become unstable.")
    print("This occurs only for oversteering vehicles, where the determinant of the system matrix becomes zero.")
    
    print("\nThe condition for an oversteering vehicle is: a * cf > b * cr")
    print("\nThe formula for critical speed is: v_crit = sqrt( (cf * cr * (a + b)^2) / (m * (a*cf - b*cr)) )")

    print("\nGiven Parameters:")
    print(f"a = {a} m")
    print(f"b = {b} m")
    print(f"cf = {cf} N/rad")
    print(f"cr = {cr} N/rad")
    print(f"m = {m} kg")

    # Step 1: Check the oversteering condition
    stability_factor = a * cf - b * cr
    
    print("\nStep 1: Check the oversteering condition (a*cf - b*cr > 0)")
    print(f"Expression: {a} * {cf} - {b} * {cr}")
    print(f"Result: {a * cf} - {b * cr} = {stability_factor}")

    if stability_factor <= 0:
        print("\nResult: The vehicle is not oversteering. A critical speed for this type of instability does not exist.")
        return

    print("The result is positive, so the vehicle is oversteering and a critical speed exists.")

    # Step 2: Calculate the terms in the formula
    wheelbase = a + b
    numerator = cf * cr * (wheelbase**2)
    denominator = m * stability_factor

    # Step 3: Calculate the critical speed
    v_crit_sq = numerator / denominator
    v_crit = math.sqrt(v_crit_sq)

    print("\nStep 2: Substitute values into the formula")
    print(f"v_crit = sqrt( ({cf} * {cr} * ({a} + {b})^2) / ({m} * ({a}*{cf} - {b}*{cr})) )")
    print(f"v_crit = sqrt( ({cf} * {cr} * ({wheelbase})**2) / ({m} * ({stability_factor})) )")
    print(f"v_crit = sqrt( {numerator:.2e} / {denominator:.2e} )")
    print(f"v_crit = sqrt( {v_crit_sq:.2f} )")
    print(f"\nFinal calculated critical speed: {v_crit:.2f} m/s")
    
    # Also print in km/h for context
    print(f"Which is equal to {v_crit * 3.6:.2f} km/h")
    
    # Return final answer in the specified format
    print(f"\n<<<{v_crit:.2f}>>>")


if __name__ == '__main__':
    calculate_critical_speed()