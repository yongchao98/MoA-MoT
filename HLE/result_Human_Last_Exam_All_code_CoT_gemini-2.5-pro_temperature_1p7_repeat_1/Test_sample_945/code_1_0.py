import math

def calculate_critical_speed(a, b, cf, cr, m):
    """
    Derives and calculates the critical speed for an oversteering vehicle.
    
    Args:
    a (float): distance from CG to front axle (m)
    b (float): distance from CG to rear axle (m)
    cf (float): front axle cornering stiffness (N/rad)
    cr (float): rear axle cornering stiffness (N/rad)
    m (float): vehicle mass (kg)
    """

    print("Derivation of Critical Speed (v_crit) for an Oversteering Vehicle\n")
    
    # Check if the vehicle is oversteering
    oversteer_term = a * cf - b * cr
    if oversteer_term <= 0:
        print("The vehicle parameters provided do not describe an oversteering vehicle.")
        print(f"(a*cf - b*cr) = {oversteer_term:.2f}, which should be positive for oversteer.")
        print("An understeering or neutral steer vehicle is stable at all speeds within this linear model.")
        return

    # Symbolic formula
    print("The formula for the square of the critical speed is:")
    print("v_crit^2 = (cf * cr * (a+b)^2) / (m * (a*cf - b*cr))\n")
    
    # Calculation
    wheelbase = a + b
    numerator = cf * cr * (wheelbase**2)
    denominator = m * oversteer_term
    
    vcrit_sq = numerator / denominator
    vcrit_ms = math.sqrt(vcrit_sq)
    vcrit_kmh = vcrit_ms * 3.6
    
    # Print the equation with substituted values
    print("Substituting the given parameter values:")
    print(f"v_crit^2 = ({cf} * {cr} * ({a} + {b})^2) / ({m} * ({a}*{cf} - {b}*{cr}))")
    print(f"v_crit^2 = ({cf * cr} * {wheelbase**2:.2f}) / ({m} * {oversteer_term})")
    print(f"v_crit^2 = {numerator:.2f} / {denominator:.2f}")
    print(f"v_crit^2 = {vcrit_sq:.2f} m^2/s^2\n")

    print("The critical speed is:")
    print(f"v_crit = sqrt({vcrit_sq:.2f}) = {vcrit_ms:.2f} m/s")
    print(f"v_crit = {vcrit_kmh:.2f} km/h")

if __name__ == '__main__':
    # --- Example Parameters for an OVERSTEERING vehicle ---
    # Note: These are example values. You can replace them with your own.
    
    # a: distance from CG to front axle (m)
    param_a = 1.2
    # b: distance from CG to rear axle (m)
    param_b = 1.6
    # cf: front axle cornering stiffness (N/rad)
    # To make it oversteer, cf is relatively high and/or cr is relatively low
    param_cf = 120000.0
    # cr: rear axle cornering stiffness (N/rad)
    param_cr = 80000.0
    # m: vehicle mass (kg)
    param_m = 1500.0
    # I: vehicle moment of inertia (kg*m^2). Note: This parameter is not needed for this specific formula but is part of the full model.

    calculate_critical_speed(param_a, param_b, param_cf, param_cr, param_m)