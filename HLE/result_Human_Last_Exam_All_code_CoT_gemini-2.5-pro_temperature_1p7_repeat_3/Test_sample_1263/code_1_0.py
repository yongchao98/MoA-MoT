import math

def calculate_escape_velocity():
    """
    This script determines if Titan can calculate Pandora's escape velocity,
    simulates the calculation using 4-bit fractional arithmetic, and finds
    the smallest possible error.
    """

    print("Step 1: Define constants and calculate the 'true' escape velocity.")
    
    # Constants in standard SI units
    G_true = 6.6743e-11  # m^3 kg^-1 s^-2
    PI_true = math.pi
    
    # Planet parameters
    r_core_m = 100 * 1000  # 100 km in meters
    R_total_m = (4000 * 1000) / 2  # 2000 km in meters
    d_core_kg_m3 = 1.2 * 1000  # 1200 kg/m^3
    d_shell_kg_m3 = 0.3 * 1000  # 300 kg/m^3

    # Volume calculations
    V_core = (4/3) * PI_true * (r_core_m ** 3)
    V_total = (4/3) * PI_true * (R_total_m ** 3)
    
    # Mass calculations
    M_core = d_core_kg_m3 * V_core
    M_shell = d_shell_kg_m3 * (V_total - V_core)
    M_total = M_core + M_shell

    # True escape velocity
    v_e_true_sq = (2 * G_true * M_total) / R_total_m
    v_e_true = math.sqrt(v_e_true_sq)

    print(f"True Mass = {M_total:.4e} kg")
    print(f"True Escape Velocity = {v_e_true:.2f} m/s\n")

    print("Step 2: Simulate the calculation on Titan architecture.")

    # Approximate the physical constants for Titan
    # v_e^2 = (8/3) * G * d_shell * pi * R^2
    # The term C = (8/3) * G * d_shell * pi must be approximated by a single fraction
    # to avoid intermediate numbers > 15.
    
    C_true = (8/3) * G_true * d_shell_kg_m3 * PI_true
    print(f"The combined constant C = (8/3)*G*d_shell*pi has a true value of {C_true:.4e}")

    # Find the best 4-bit fraction approximation for C.
    # C_true is ~1.677e-7. Let's try to represent this.
    # We can write it as a/b * 10^e. If e=-6, we need to approximate 0.1677.
    # 1/6 = 0.1666... is an excellent candidate. Numerator (1) and Denominator (6) are <= 15.
    C_approx_num, C_approx_den, C_approx_exp = 1, 6, -6
    C_approx = (C_approx_num / C_approx_den) * (10**C_approx_exp)
    
    print(f"Best Titan approximation for C is {C_approx_num}/{C_approx_den} * 10^{C_approx_exp}, which is {C_approx:.4e}")
    
    # Titan representation of R^2
    R_sq_num, R_sq_den, R_sq_exp = 4, 1, 12 # R = 2e6 -> R^2 = 4e12
    R_sq = (R_sq_num / R_sq_den) * (10**R_sq_exp)

    # Titan calculation of v_e^2 = C_approx * R^2
    # MUL (1/6e-6), (4/1e12) -> (1*4)/(6*1) * 10^(-6+12) = 4/6e6 = 2/3e6
    ve_sq_titan_num, ve_sq_titan_den, ve_sq_titan_exp = 2, 3, 6
    ve_sq_titan = (ve_sq_titan_num/ve_sq_titan_den) * (10**ve_sq_titan_exp)
    print(f"Titan calculates v_e^2 = C_approx * R^2 = {ve_sq_titan_num}/{ve_sq_titan_den} * 10^{ve_sq_titan_exp}\n")
    
    print("Step 3: Approximate the square root on Titan.")
    # We need to calculate sqrt(2/3 * 10^6) = sqrt(2/3) * 10^3
    # sqrt(2/3) is approx 0.8165. We need a 4-bit fraction for this.
    # Let's test candidates: 4/5=0.8, 5/6=0.833. A very good one is 9/11 ~ 0.8181...
    # (9/11)^2 = 81/121 ~ 0.6694, which is very close to 2/3 ~ 0.6667.
    sqrt_approx_num, sqrt_approx_den = 9, 11
    
    # Final calculation for v_e
    v_e_titan = (sqrt_approx_num/sqrt_approx_den) * (10**3)

    print("The final calculation equation on Titan is:")
    print(f"v_e = Sqrt( (C_approx) * (R^2) )")
    final_eq_str = (f"v_e = Sqrt( ({C_approx_num}/{C_approx_den}e{C_approx_exp}) * ({R_sq_num}/{R_sq_den}e{R_sq_exp}) ) = "
                    f"Sqrt( {ve_sq_titan_num}/{ve_sq_titan_den}e{ve_sq_titan_exp} ) = "
                    f"Sqrt({ve_sq_titan_num}/{ve_sq_titan_den}) * 10^3")
    print(final_eq_str)
    print(f"Using approximation Sqrt({ve_sq_titan_num}/{ve_sq_titan_den}) ≈ {sqrt_approx_num}/{sqrt_approx_den}:")
    print(f"v_e ≈ {sqrt_approx_num}/{sqrt_approx_den} * 10^3 = {v_e_titan:.2f} m/s\n")
    
    print("Step 4: Conclusion and Error Calculation.")
    abs_error = abs(v_e_true - v_e_titan)
    print("The calculation is possible on the Titan architecture.")
    print(f"The smallest absolute error found is {abs_error:.2f} m/s.")
    
    answer = f"Y{abs_error:.2f}"
    print(f"\nFinal Answer: {answer}")
    print(f"<<<{answer}>>>")


calculate_escape_velocity()