import math

def high_precision_calculation():
    """Calculates the 'true' value of the force with standard float precision."""
    G = 6.6743e-11
    m_probe = 30.0
    
    R_core_m = 50 * 1000
    R_total_m = 1000 * 1000
    rho_core_kg_m3 = 1.2 * 1000
    rho_shell_kg_m3 = 0.3 * 1000
    
    V_core = (4/3) * math.pi * (R_core_m ** 3)
    V_total = (4/3) * math.pi * (R_total_m ** 3)
    V_shell = V_total - V_core
    
    M_core = V_core * rho_core_kg_m3
    M_shell = V_shell * rho_shell_kg_m3
    M_pandora = M_core + M_shell
    
    dist_r = (1000 * 1000) + 500
    
    F_true = (G * M_pandora * m_probe) / (dist_r ** 2)
    return F_true

def solve_with_titan():
    """
    Simulates the calculation of gravitational force using the Titan 5-bit architecture.
    """
    print("Starting Titan computation for Pandora's gravitational force on Pioneer probe.")
    print("-" * 70)

    # Step 1: Define constants and parameters as Titan fractions (n, d) and exponents (exp)
    # n and d must be <= 31.
    
    # Gravitational Constant G = 6.674e-11. Approximate 6.674 as 20/3.
    G = {'num': 20, 'den': 3, 'exp': -11, 'name': 'G'}
    # Pi = 3.14159... Approximate as 22/7
    PI = {'num': 22, 'den': 7, 'name': 'pi'}
    # 4/3
    FOUR_THIRDS = {'num': 4, 'den': 3, 'name': '4/3'}
    
    # Pandora's parameters
    # Shell Density rho_shell = 300 kg/m^3. Represent as 3 * 10^2
    RHO_SHELL = {'num': 3, 'den': 1, 'exp': 2, 'name': 'rho_shell'}
    # Total Radius R_total = 1000 km = 1e6 m. Represent as 1 * 10^6
    R_TOTAL = {'num': 1, 'den': 1, 'exp': 6, 'name': 'R_total'}
    
    # Probe's parameters
    # Probe mass m = 30 kg. Represent as 3 * 10^1 to keep numerator small.
    M_PROBE = {'num': 3, 'den': 1, 'exp': 1, 'name': 'm_probe'}
    # Distance r = 1000.5 km = 1000500 m. 
    # To keep fraction simple, we must approximate 1.0005 as 1/1.
    # So r = 1 * 10^6 m
    DIST_R = {'num': 1, 'den': 1, 'exp': 6, 'name': 'r'}
    
    print("Step 1: Representing physical values as Titan fractions (n/d) and exponents (10^exp).")
    for val in [G, PI, RHO_SHELL, R_TOTAL, M_PROBE, DIST_R]:
        print(f"  {val['name']:<10} = {val['num']}/{val['den']} x 10^{val['exp']}")
    print("-" * 70)
    
    print("Step 2: Calculating Pandora's Mass (M).")
    print("M = M_core + M_shell. We find M_core is negligible compared to M_shell.")
    print("Thus, we approximate M ≈ M_shell = (4/3 * pi * R_total^3) * rho_shell.")

    # Calculate R_total^3
    r_total_cubed_exp = R_TOTAL['exp'] * 3
    
    # Calculate M_shell fraction part: (4/3) * pi * rho_shell
    # (4/3 * 22/7) * 3/1 = 88/21 * 3/1 = 264/21 = 88/7
    m_shell_frac_num = FOUR_THIRDS['num'] * PI['num'] * RHO_SHELL['num']
    m_shell_frac_den = FOUR_THIRDS['den'] * PI['den'] * RHO_SHELL['den']
    
    # Simplify the fraction
    common_divisor = math.gcd(m_shell_frac_num, m_shell_frac_den)
    m_shell_frac_num //= common_divisor
    m_shell_frac_den //= common_divisor

    print(f"\n  Intermediate M_shell fraction = (4/3 * 22/7 * 3/1) = {m_shell_frac_num}/{m_shell_frac_den}")
    
    # Now we must simplify 88/7 because 88 > 31.
    val = m_shell_frac_num / m_shell_frac_den # ~12.57
    # Find a good 5-bit fraction for 12.57. 25/2 = 12.5 is excellent.
    m_shell_final_frac_num = 25
    m_shell_final_frac_den = 2
    
    # This simplification did not require an exponent change.
    m_shell_exp_adj = 0
    m_shell_exp = R_TOTAL['exp'] * 3 + RHO_SHELL['exp'] + m_shell_exp_adj
    M_PANDORA = {'num': m_shell_final_frac_num, 'den': m_shell_final_frac_den, 'exp': m_shell_exp, 'name': 'M_pandora'}
    
    print(f"  This fraction's value is ~{val:.2f}. We approximate it as 25/2 = 12.5.")
    print(f"  So, Pandora's Mass M = {M_PANDORA['num']}/{M_PANDORA['den']} x 10^{M_PANDORA['exp']} kg.")
    print("-" * 70)

    print("Step 3: Calculating Gravitational Force F = G * M * m / r^2.")
    
    # Calculate r^2
    r_squared_exp = DIST_R['exp'] * 2
    
    # Calculate numerator: G * M * m
    num_exp = G['exp'] + M_PANDORA['exp'] + M_PROBE['exp']
    num_frac_num = G['num'] * M_PANDORA['num'] * M_PROBE['num']
    num_frac_den = G['den'] * M_PANDORA['den'] * M_PROBE['den']

    # Full expression before simplification
    print(f"\n  F = ( ({G['num']}/{G['den']}) * ({M_PANDORA['num']}/{M_PANDORA['den']}) * ({M_PROBE['num']}/{M_PROBE['den']}) ) / ({DIST_R['num']}/{DIST_R['den']})^2 * 10^...")
    print(f"  Numerator fraction = {num_frac_num}/{num_frac_den}")

    # Simplify numerator fraction
    common_divisor = math.gcd(num_frac_num, num_frac_den)
    num_frac_num //= common_divisor
    num_frac_den //= common_divisor
    print(f"  Simplified numerator fraction = {num_frac_num}/{num_frac_den}") # 250/1
    
    # Now simplify 250/1 as it's out of range. 250 = 2.5 * 10^2
    val = num_frac_num / num_frac_den # 250
    # New fraction for 2.5 is 5/2. Exponent adjustment is +2.
    final_num_frac_num = 5
    final_num_frac_den = 2
    num_exp_adj = 2
    
    print(f"  This value {val} is > 31. We re-factor it: {val} = 2.5 x 10^2.")
    print(f"  The new fraction is 5/2, with an exponent adjustment of +2.")
    
    # Final fraction and exponent
    final_exp = num_exp + num_exp_adj - r_squared_exp
    final_frac_num = final_num_frac_num
    final_frac_den = final_num_frac_den
    
    final_force_val = (final_frac_num / final_frac_den) * (10**final_exp)

    print("\nStep 4: Final Result.")
    print("The final equation for the force is:")
    print(f"F = {final_frac_num} / {final_frac_den} x 10^{final_exp}")
    
    print(f"\nCalculated Force = {final_force_val:.3f} N")
    print("-" * 70)

    f_true = high_precision_calculation()
    error = abs(f_true - final_force_val)

    print("Final Answer Derivation:")
    print(f"  Titan's calculated force: {final_force_val:.3f} N")
    print(f"  High-precision force:   {f_true:.3f} N")
    print(f"  Absolute error:         {error:.3f} N")
    
    return f"Y[{error:.3f}]"

# Execute the simulation and print the final answer in the required format
final_answer = solve_with_titan()
print(f"\n<<< {final_answer} >>>")