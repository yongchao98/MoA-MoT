import math

def calculate_critical_speed(a, b, c_f, c_r, m):
    """
    Calculates the critical speed for an oversteering vehicle using the linear single-track model.

    Args:
        a (float): Distance from CG to front axle (m).
        b (float): Distance from CG to rear axle (m).
        c_f (float): Cornering stiffness of the front axle (N/rad).
        c_r (float): Cornering stiffness of the rear axle (N/rad).
        m (float): Vehicle mass (kg).
    """
    print("Deriving critical speed for a single-track vehicle model.")
    print("-" * 50)
    print("Parameters:")
    print(f"a (CG to front axle) = {a} m")
    print(f"b (CG to rear axle) = {b} m")
    print(f"c_f (front cornering stiffness) = {c_f} N/rad")
    print(f"c_r (rear cornering stiffness) = {c_r} N/rad")
    print(f"m (vehicle mass) = {m} kg")
    print("-" * 50)

    # Check for oversteer condition
    oversteer_term = a * c_f - b * c_r
    if oversteer_term <= 0:
        print("Vehicle is not oversteering (ac_f - bc_r <= 0).")
        print("Based on this model, a real critical speed for this type of instability does not exist.")
        return

    # Calculate critical speed
    L = a + b
    numerator = (L**2) * c_f * c_r
    denominator = m * oversteer_term

    v_squared = numerator / denominator
    v_crit = math.sqrt(v_squared)

    # Print the equation step-by-step
    print("The critical speed (v_crit) is calculated using the formula:")
    print("v_crit = sqrt( ((a+b)^2 * c_f * c_r) / (m * (a*c_f - b*c_r)) )\n")

    print("Substituting the values:")
    print(f"Numerator = (({a} + {b})^2 * {c_f} * {c_r}) = {numerator}")
    print(f"Denominator = ({m} * ({a} * {c_f} - {b} * {c_r})) = {denominator}")
    print("\nCalculation:")
    print(f"v_crit = sqrt( {numerator} / {denominator} )")
    print(f"v_crit = sqrt( {v_squared:.2f} )")
    
    print("-" * 50)
    print(f"Final Result: Critical Speed (v_crit) = {v_crit:.2f} m/s")
    print("-" * 50)
    
    return v_crit

if __name__ == '__main__':
    # --- Vehicle Parameters ---
    # These parameters describe an oversteering vehicle (a*c_f - b*c_r > 0)
    # The condition for oversteer is a/b > c_r/c_f
    a_param = 1.8  # m, distance from CG to front axle
    b_param = 0.9  # m, distance from CG to rear axle
    cf_param = 80000.0  # N/rad, cornering stiffness of front axle
    cr_param = 70000.0  # N/rad, cornering stiffness of rear axle
    m_param = 1500.0  # kg, vehicle mass
    # I_param is not needed for this calculation but is part of the full model.

    calculate_critical_speed(a_param, b_param, cf_param, cr_param, m_param)