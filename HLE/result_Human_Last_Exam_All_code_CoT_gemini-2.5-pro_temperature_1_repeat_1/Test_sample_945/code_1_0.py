import math

def derive_critical_speed():
    """
    Calculates and prints the critical speed for an oversteering vehicle.
    
    The function uses example parameters for an oversteering vehicle and prints the
    derivation formula, the formula with values, and the final result.
    """
    # Parameters of the linear single-track model
    # a: distance from CG to front axle [m]
    # b: distance from CG to rear axle [m]
    # c_f: cornering stiffness of the front axle [N/rad]
    # c_r: cornering stiffness of the rear axle [N/rad]
    # m: vehicle mass [kg]
    
    # Example parameters for an OVERSTEERING vehicle (where a*c_f > b*c_r)
    a = 1.5   # m
    b = 1.0   # m
    c_f = 80000 # N/rad
    c_r = 60000 # N/rad
    m = 1500  # kg

    # The moment of inertia 'I' is not needed for the stability boundary calculation.

    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle\n")
    
    # Check for the oversteering condition
    oversteer_term = a * c_f - b * c_r
    if oversteer_term <= 0:
        print("The provided parameters do not describe an oversteering vehicle.")
        print(f"Condition (a*c_f - b*c_r > 0) is not met: {oversteer_term:.2f} <= 0")
        return

    # Print the symbolic formula
    print("Formula:")
    print("v_crit = sqrt( ((a+b)^2 * c_f * c_r) / (m * (a*c_f - b*c_r)) )\n")

    # Print the formula with the numerical values plugged in
    print("Calculation with example values:")
    L = a + b
    numerator_str = f"(({a} + {b})^2 * {c_f} * {c_r})"
    denominator_str = f"({m} * ({a}*{c_f} - {b}*{c_r}))"
    print(f"v_crit = sqrt( {numerator_str} / {denominator_str} )\n")

    # Perform the calculation
    numerator_val = (L**2) * c_f * c_r
    denominator_val = m * oversteer_term
    
    v_crit_sq = numerator_val / denominator_val
    v_crit_ms = math.sqrt(v_crit_sq) # result in meters per second
    v_crit_kmh = v_crit_ms * 3.6     # convert to kilometers per hour

    # Print the final result
    print("Result:")
    print(f"The critical speed is {v_crit_ms:.2f} m/s, which is equivalent to {v_crit_kmh:.2f} km/h.")
    print("Above this speed, the vehicle's linear lateral dynamics become unstable.")


if __name__ == '__main__':
    derive_critical_speed()
    # The final numerical answer for the example parameters is ~18.26 m/s
    # We output this via the print statements in the function.
    # To satisfy the format requirement, let's output the m/s value.
    final_answer = math.sqrt( ((1.5+1.0)**2 * 80000 * 60000) / (1500 * (1.5*80000 - 1.0*60000)) )
    # print(f'<<<{final_answer:.2f}>>>') #This would be printed if needed by the system
    
# The final answer is derived and printed within the function.
# The user can just execute the script.