import math

def find_best_fraction(target_val, max_int=31):
    """Finds the best fractional approximation a/b for a target value."""
    best_frac = (0, 1)
    min_error = float('inf')
    for b in range(1, max_int + 1):
        # Find the closest integer numerator 'a' for the given 'b'
        a_ideal = target_val * b
        a_low = math.floor(a_ideal)
        a_high = math.ceil(a_ideal)

        # Check the lower candidate
        if a_low >= 0 and a_low <= max_int:
            error = abs(target_val - a_low / b)
            if error < min_error:
                min_error = error
                best_frac = (a_low, b)

        # Check the higher candidate
        if a_high >= 0 and a_high <= max_int and a_high != a_low:
            error = abs(target_val - a_high / b)
            if error < min_error:
                min_error = error
                best_frac = (a_high, b)
    return best_frac

def main():
    """
    Solves the Pandora landing problem using Titan 5-bit architecture simulation.
    """
    # Step 1: Define constants and calculate "true" values
    G_true = 6.674e-11
    R_core_true = 1e5  # m
    R_pandora_true = 2e6 # m
    rho_core_true = 1200 # kg/m^3
    rho_shell_true = 300 # kg/m^3
    
    # M = (4/3) * pi * [r_c^3 * rho_c + (R^3 - r_c^3) * rho_s]
    M_pandora_true = (4/3) * math.pi * (R_core_true**3 * rho_core_true + (R_pandora_true**3 - R_core_true**3) * rho_shell_true)
    
    m_probe = 50.0  # kg
    v_initial = 300.0 # m/s
    d_landing = 5000.0 # m
    r_probe_true = R_pandora_true + d_landing # Probe altitude is negligible compared to radius, use R_pandora
    
    g_true = G_true * M_pandora_true / R_pandora_true**2
    a_decel_true = v_initial**2 / (2 * d_landing)
    
    F_g_true = m_probe * g_true
    F_decel_true = m_probe * a_decel_true
    F_rocket_true = F_g_true + F_decel_true

    print("--- Titan Architecture Calculation ---")
    print("Objective: Calculate F_rocket = F_gravity + F_deceleration\n")

    # Step 2: Represent initial values in Titan format
    # G ≈ 20/3 * 10^-11 (error ~0.1%)
    G_titan = ((20, 3), -11)
    # M ≈ 1 * 10^22 kg. Based on M_true ≈ 1.0056e22 kg.
    M_titan = ((1, 1), 22)
    # R^2 = (2e6)^2 = 4e12
    R2_titan = ((4, 1), 12)
    # m = 50 kg = 5 * 10^1
    m_titan = ((5, 1), 1)
    # v = 300 m/s = 3 * 10^2
    v_titan = ((3, 1), 2)
    # d = 5000 m = 5 * 10^3
    d_titan = ((5, 1), 3)

    # Step 3: Calculate F_gravity
    print("1. Calculate Gravitational Force (F_g = m * g)")
    # g = G*M/R^2 = ((20/3)e-11 * 1e22) / 4e12
    # Numerator: (20/3) * 1 = 20/3. Exp: -11+22=11. -> (20/3)*10^11
    # Denominator: 4*10^12
    # Division: (20/3) / 4 = 20/12 = 5/3. Exp: 11-12=-1
    g_titan_frac = (5, 3) # represents 1.666...
    g_titan_exp = -1
    print(f"   - Calculated g = ({g_titan_frac[0]}/{g_titan_frac[1]}) * 10^{g_titan_exp} m/s^2")

    # F_g = m * g = (5*10^1) * (5/3 * 10^-1)
    # Frac: 5 * 5/3 = 25/3. Exp: 1-1=0.
    # 25 and 3 are valid 5-bit integers.
    Fg_titan_frac = (25, 3)
    Fg_titan_exp = 0
    print(f"   - F_g = m * g = ({Fg_titan_frac[0]}/{Fg_titan_frac[1]}) * 10^{Fg_titan_exp} N\n")

    # Step 4: Calculate F_deceleration
    print("2. Calculate Deceleration Force (F_decel = m * a)")
    # a = v^2 / (2d)
    # v^2 = (3*10^2)^2 = 9 * 10^4
    v2_titan_frac = (9, 1)
    v2_titan_exp = 4
    # 2d = 2 * (5*10^3) = 10 * 10^3 = 1 * 10^4
    two_d_titan_frac = (1, 1)
    two_d_titan_exp = 4
    # a = (9*10^4) / (1*10^4) = 9
    a_decel_titan_frac = (9, 1)
    a_decel_titan_exp = 0
    print(f"   - Calculated a = ({a_decel_titan_frac[0]}/{a_decel_titan_frac[1]}) * 10^{a_decel_titan_exp} m/s^2")

    # F_decel = m * a = (5*10^1) * 9
    # The multiplication (5/1)*(9/1) results in 45/1. Numerator 45 > 31, so this is invalid.
    # We must calculate the value first and then represent it.
    # F_decel = 50 * 9 = 450 N.
    # Represent 450 in Titan format: 4.5 * 10^2 = (9/2) * 10^2
    Fdecel_titan_frac = (9, 2)
    Fdecel_titan_exp = 2
    print(f"   - F_decel = 450 N, represented as ({Fdecel_titan_frac[0]}/{Fdecel_titan_frac[1]}) * 10^{Fdecel_titan_exp} N\n")
    
    # Step 5: Add the two forces
    print("3. Add Forces (F_rocket = F_g + F_decel)")
    # F_g = 25/3 = 8.333... N
    # F_decel = (9/2) * 10^2 = 450 N
    # To add, exponents must match. Convert F_g to have exponent 2.
    # F_g = 8.333... = 0.08333... * 10^2
    # We need to approximate 0.08333... (which is 1/12)
    Fg_adj_frac = (1, 12)
    Fg_adj_exp = 2
    print(f"   - To add, convert F_g to common exponent: {Fg_titan_frac[0]}/{Fg_titan_frac[1]} N -> ({Fg_adj_frac[0]}/{Fg_adj_frac[1]}) * 10^{Fg_adj_exp} N")

    # Add fractional parts: 1/12 + 9/2
    # 1/12 + 54/12 = 55/12
    sum_num = Fg_adj_frac[0] * Fdecel_titan_frac[1] + Fdecel_titan_frac[0] * Fg_adj_frac[1]
    sum_den = Fg_adj_frac[1] * Fdecel_titan_frac[1]
    print(f"   - Sum of fractions: {Fg_adj_frac[0]}/{Fg_adj_frac[1]} + {Fdecel_titan_frac[0]}/{Fdecel_titan_frac[1]} = {sum_num}/{sum_den}")

    # Numerator 55 exceeds 31. Must approximate.
    target_value = sum_num / sum_den
    print(f"   - Numerator {sum_num} > 31. Approximating {sum_num}/{sum_den} ({target_value:.4f})...")
    final_frac = find_best_fraction(target_value)
    
    F_rocket_titan = (final_frac[0]/final_frac[1]) * (10**Fdecel_titan_exp)
    print(f"   - Best approximation with 5-bit integers: {final_frac[0]}/{final_frac[1]} ({final_frac[0]/final_frac[1]:.4f})\n")

    print("--- Final Result ---")
    print("The final Titan calculation for the force is:")
    print(f"F_rocket ≈ ({Fg_adj_frac[0]}/{Fg_adj_frac[1]})*10^{Fg_adj_exp} + ({Fdecel_titan_frac[0]}/{Fdecel_titan_frac[1]})*10^{Fdecel_titan_exp} ≈ ({final_frac[0]}/{final_frac[1]})*10^{Fdecel_titan_exp}")
    
    print(f"\nThe numbers in the final step of the equation are: {Fg_adj_frac[0]}, {Fg_adj_frac[1]}, {Fdecel_titan_frac[0]}, {Fdecel_titan_frac[1]}, {final_frac[0]}, {final_frac[1]}, and the exponent {Fdecel_titan_exp}.")

    print(f"\nCalculated Titan Force = {F_rocket_titan:.3f} N")
    print(f"Baseline 'True' Force = {F_rocket_true:.3f} N")
    
    abs_error = abs(F_rocket_titan - F_rocket_true)
    print(f"Smallest Absolute Error (e) = {abs_error:.3f}")

    print("\nAnswer: Yes, the force can be calculated.")
    print(f"<<<Y[{abs_error:.3f}]>>>")

if __name__ == "__main__":
    main()