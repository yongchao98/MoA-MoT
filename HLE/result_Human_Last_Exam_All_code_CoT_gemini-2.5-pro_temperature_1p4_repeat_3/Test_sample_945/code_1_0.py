import math

def calculate_critical_speed():
    """
    Derives and calculates the critical speed for an oversteering vehicle.
    """
    # Parameters of the linear single-track model (example values)
    a = 1.5      # distance CG to front axle (m)
    b = 1.0      # distance CG to rear axle (m)
    c_f = 90000  # cornering stiffness of front axle (N/rad)
    c_r = 100000 # cornering stiffness of rear axle (N/rad)
    m = 1500     # vehicle mass (kg)
    I = 2500     # vehicle moment of inertia (kg*m^2)

    # The derived formula for the square of the critical speed is:
    # v_crit^2 = (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r))
    # Stability is lost when the vehicle is oversteering, which corresponds to the
    # condition where the denominator (a * c_f - b * c_r) is positive.

    print("--- Derivation of Critical Speed ---")
    print("The derived formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( (c_f * c_r * (a + b)**2) / (m * (a * c_f - b * c_r)) )")
    print("-" * 35)

    # Check for the oversteering condition
    oversteer_term = a * c_f - b * c_r

    if oversteer_term > 0:
        print("Vehicle parameters satisfy the oversteering condition (a*c_f > b*c_r).")
        print("A critical speed exists.\n")

        # Calculate numerator and denominator
        numerator = c_f * c_r * (a + b)**2
        denominator = m * oversteer_term
        
        # Calculate critical speed
        v_crit_sq = numerator / denominator
        v_crit = math.sqrt(v_crit_sq)

        # Output the calculation with numbers
        print("--- Calculation with Provided Parameters ---")
        print(f"a = {a} m, b = {b} m, c_f = {c_f} N/rad, c_r = {c_r} N/rad, m = {m} kg, I = {I} kg*m^2")
        print("\nSubstituting the values into the formula:")
        
        # Show each number in the equation
        print(f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})**2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )")
        print(f"v_crit = sqrt( ({c_f * c_r} * {(a + b)**2}) / ({m} * {oversteer_term}) )")
        print(f"v_crit = sqrt( {numerator:.2e} / {denominator:.2e} )")
        print(f"v_crit = sqrt( {v_crit_sq:.2f} m^2/s^2 )")
        print(f"\nCritical Speed = {v_crit:.2f} m/s")
        print(f"Critical Speed = {v_crit * 3.6:.2f} km/h")

    else:
        print("Vehicle is not oversteering (a*c_f <= b*c_r).")
        print("The linear model is stable at all speeds.")

# Run the calculation
calculate_critical_speed()