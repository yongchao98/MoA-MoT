import math

def calculate_pandora_escape_velocity():
    """
    This function calculates the escape velocity of Pandora using both
    high-precision floating-point arithmetic for the "true" value and simulates
    the best possible approximation on the Titan architecture.
    """
    # 1. High-precision constants and planet parameters
    G_TRUE = 6.6743e-11  # m^3 kg^-1 s^-2
    PI_TRUE = math.pi

    # Core parameters
    r_c = 100 * 1000  # 100 km to meters
    rho_c = 1.2 * 1000 # 1.2 metric tons/m^3 to kg/m^3

    # Total parameters
    R = 4000 * 1000 / 2 # 4000 km diameter to radius in meters
    rho_s = 0.3 * 1000  # 0.3 metric tons/m^3 to kg/m^3

    # 2. Calculate true mass and escape velocity
    v_c = (4/3) * PI_TRUE * (r_c ** 3)
    m_c = rho_c * v_c

    v_total = (4/3) * PI_TRUE * (R ** 3)
    v_s = v_total - v_c
    m_s = rho_s * v_s

    m_total = m_c + m_s
    ve_true = math.sqrt(2 * G_TRUE * m_total / R)

    print("Step 1: Calculating the true escape velocity of Pandora.")
    print(f"  - Pandora's Total Mass M = {m_total:.4e} kg")
    print(f"  - Pandora's Radius R = {R:.4e} m")
    print(f"  - Gravitational Constant G = {G_TRUE:.4e} m^3 kg^-1 s^-2")
    print(f"  - The true escape velocity is calculated as sqrt(2*G*M/R).")
    print(f"  - True Escape Velocity (ve_true) = {ve_true:.4f} m/s\n")

    # 3. Find the best Titan approximation for the final result
    print("Step 2: Simulating Titan's RED instruction to find the best fractional approximation.")
    print("  - We search for the best fraction (n/d) and exponent (k) to represent the value.")
    
    min_error = float('inf')
    best_approx_val = 0
    best_approx_details = {}

    # Titan's integers are 4-bit (0-15)
    INT_MAX = 15
    
    # Search over a reasonable range of exponents for scientific notation
    for k in range(5): # searching for 10^0, 10^1, ..., 10^4
        target = ve_true / (10**k)
        
        for d in range(1, INT_MAX + 1):
            for n in range(0, INT_MAX + 1):
                # Avoid division by zero if n=d=0, though loop starts at d=1
                current_val_base = n / d
                error = abs(current_val_base - target)
                
                if error < min_error:
                    min_error = error
                    best_approx_val = current_val_base * (10**k)
                    best_approx_details = {
                        'n': n,
                        'd': d,
                        'k': k,
                        'val': best_approx_val
                    }

    final_abs_error = abs(ve_true - best_approx_val)

    print(f"  - The best approximation found is ({best_approx_details['n']}/{best_approx_details['d']}) * 10^{best_approx_details['k']}.")
    print(f"  - Titan-approximated velocity = {best_approx_details['val']:.4f} m/s\n")
    
    print("Step 3: Calculating the smallest absolute error.")
    print(f"  - Error = | True Value - Titan Value |")
    print(f"  - Error = | {ve_true:.4f} - {best_approx_details['val']:.4f} |")
    print(f"  - Final Absolute Error = {final_abs_error:.4f} m/s")
    print(f"  - Rounded to 0.01 m/s, the error is {final_abs_error:.2f} m/s.\n")

    # The final answer format required
    final_answer = f"Y[{final_abs_error:.2f}]"
    print(f"The calculation is possible. The answer in the required format is: <<< {final_answer} >>>")


calculate_pandora_escape_velocity()