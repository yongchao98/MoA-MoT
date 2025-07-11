import math

def calculate_critical_speed():
    """
    Calculates the critical speed for an oversteering vehicle based on the linear single-track model.
    """
    # --- Vehicle Parameters (replace with your values) ---
    # a: distance from CG to front axle (m)
    a = 1.2
    # b: distance from CG to rear axle (m)
    b = 1.6
    # c_f: cornering stiffness of the front axle (N/rad)
    c_f = 80000.0
    # c_r: cornering stiffness of the rear axle (N/rad)
    c_r = 120000.0
    # m: vehicle mass (kg)
    m = 1500.0
    
    # Check for oversteering condition
    oversteer_check = b * c_r - a * c_f
    if oversteer_check <= 0:
        print("The provided parameters do not describe an oversteering vehicle.")
        print(f"The term (b*c_r - a*c_f) must be positive, but it is: {oversteer_check:.2f}")
        return

    # Calculate the numerator of the fraction
    # numerator = (a + b)^2 * c_f * c_r
    numerator = (a + b)**2 * c_f * c_r

    # Calculate the denominator of the fraction
    # denominator = m * (b*c_r - a*c_f)
    denominator = m * oversteer_check

    # Calculate the critical speed squared
    v_crit_sq = numerator / denominator

    # Calculate the critical speed (m/s)
    v_crit_ms = math.sqrt(v_crit_sq)
    
    # Convert to km/h for reference
    v_crit_kmh = v_crit_ms * 3.6

    # Print the equation with the values substituted
    print("Derivation of critical speed (v_crit):")
    print(f"v_crit = sqrt( (({a} + {b})^2 * {c_f} * {c_r}) / ({m} * ({b} * {c_r} - {a} * {c_f})) )")
    print(f"v_crit = sqrt( (({a + b})**2 * {c_f} * {c_r}) / ({m} * ({oversteer_check})) )")
    print(f"v_crit = sqrt( ({numerator}) / ({denominator}) )")
    print(f"v_crit = sqrt({v_crit_sq:.4f})")
    
    # Print the final result
    print("\n--- Result ---")
    print(f"The critical speed is: {v_crit_ms:.2f} m/s")
    print(f"Which is equivalent to: {v_crit_kmh:.2f} km/h")

# Execute the calculation
calculate_critical_speed()