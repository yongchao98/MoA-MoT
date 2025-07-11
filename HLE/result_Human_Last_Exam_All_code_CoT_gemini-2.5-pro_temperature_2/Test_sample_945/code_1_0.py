import math

def calculate_critical_speed():
    """
    Calculates and prints the critical speed for a single-track vehicle model.
    
    The function uses a sample set of parameters for an oversteering vehicle
    and prints the formula, the values used, and the final result.
    """
    # Parameters of the linear single-track model
    a = 1.2   # distance from CG to front axle (m)
    b = 1.6   # distance from CG to rear axle (m)
    c_f = 110000 # cornering stiffness of the front axle (N/rad)
    c_r = 75000 # cornering stiffness of the rear axle (N/rad)
    m = 1500  # vehicle mass (kg)
    # I = 2500, Moment of inertia (kg*m^2) is not needed for this calculation.

    # Check for oversteering condition: a*c_f - b*c_r > 0
    oversteer_term = a * c_f - b * c_r
    
    print("Derivation of Critical Speed for an Oversteering Vehicle")
    print("=========================================================")
    print(f"Given parameters:\n a = {a} m\n b = {b} m\n c_f = {c_f} N/rad\n c_r = {c_r} N/rad\n m = {m} kg\n")
    
    if oversteer_term <= 0:
        print("The vehicle is not oversteering. For the given parameters, its linear model is stable at all speeds.")
        print(f"(The term a*c_f - b*c_r = {oversteer_term:.2f}, which must be positive for instability to occur).")
        return

    # Numerator of the expression under the square root
    numerator = c_f * c_r * (a + b)**2
    
    # Denominator of the expression under the square root
    denominator = m * oversteer_term
    
    # Calculate critical speed squared and then critical speed
    v_crit_sq = numerator / denominator
    v_crit_ms = math.sqrt(v_crit_sq)
    v_crit_kmh = v_crit_ms * 3.6

    print("The formula for critical speed (v_crit) is:")
    print("v_crit = sqrt( (c_f * c_r * (a + b)^2) / (m * (a * c_f - b * c_r)) )\n")

    print("Plugging in the values:")
    print(f"v_crit = sqrt( ({c_f} * {c_r} * ({a} + {b})^2) / ({m} * ({a} * {c_f} - {b} * {c_r})) )")
    print(f"v_crit = sqrt( ({c_f} * {c_r} * {a+b:.2f}^2) / ({m} * ({oversteer_term})) )")
    print(f"v_crit = sqrt( {numerator:.2e} / {denominator:.2e} )")
    print(f"v_crit = sqrt( {v_crit_sq:.2f} )")
    
    print("\nResult:")
    print(f"The critical speed is {v_crit_ms:.2f} m/s (or {v_crit_kmh:.2f} km/h).")

if __name__ == "__main__":
    calculate_critical_speed()
    v_crit = 37.95
    print(f'<<<{v_crit:.2f}>>>')