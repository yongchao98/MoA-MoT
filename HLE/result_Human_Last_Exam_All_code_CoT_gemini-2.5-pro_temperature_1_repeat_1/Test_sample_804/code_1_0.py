import math

def simplify_fractions(numerators, denominators):
    """
    Simplifies a product of fractions by iteratively finding and canceling common factors.
    """
    # Make copies to modify
    nums = list(numerators)
    dens = list(denominators)
    
    changed = True
    while changed:
        changed = False
        for i, n in enumerate(nums):
            for j, d in enumerate(dens):
                common_divisor = math.gcd(n, d)
                if common_divisor > 1:
                    nums[i] //= common_divisor
                    dens[j] //= common_divisor
                    changed = True
            if changed:
                break
    
    # Calculate the final numerator and denominator
    final_num = 1
    for n in nums:
        final_num *= n
        
    final_den = 1
    for d in dens:
        final_den *= d
        
    return final_num, final_den, nums, dens

def main():
    """
    Performs the gravitational force calculation based on Titan's architecture.
    """
    # 1. Define constants and variables in Titan's fractional format (value, exponent)
    # G = 6.674e-11 N m^2/kg^2. Approximate as 20/3 = 6.667
    G_frac = (20, 3)
    G_exp = -11
    
    # m_p = 30 kg
    m_p_frac = (30, 1)
    m_p_exp = 0
    
    # 4/3
    four_thirds_frac = (4, 3)
    four_thirds_exp = 0

    # pi = 3.14159... Approximate as 25/8 = 3.125 to allow for calculation
    pi_frac = (25, 8)
    pi_exp = 0
    
    # r_s = 1000 km = 1e6 m
    r_s_frac = (1, 1)
    r_s_exp = 6
    
    # rho_s = 0.3 tons/m^3 = 300 kg/m^3
    rho_s_frac = (3, 1)
    rho_s_exp = 2
    
    # 2. Combine fractions and exponents from the simplified formula F ≈ G * m_p * (4/3) * pi * r_s * rho_s
    
    numerators = [G_frac[0], m_p_frac[0], four_thirds_frac[0], pi_frac[0], r_s_frac[0], rho_s_frac[0]]
    denominators = [G_frac[1], m_p_frac[1], four_thirds_frac[1], pi_frac[1], r_s_frac[1], rho_s_frac[1]]
    
    total_exponent = G_exp + m_p_exp + four_thirds_exp + pi_exp + r_s_exp + rho_s_exp

    print("Calculating the force F ≈ G * m_p * (4/3) * π * r_s * ρ_s")
    print("Equation with Titan's fractional approximations:")
    print(f"F ≈ ({G_frac[0]}/{G_frac[1]}) * ({m_p_frac[0]}/{m_p_frac[1]}) * ({four_thirds_frac[0]}/{four_thirds_frac[1]}) * ({pi_frac[0]}/{pi_frac[1]}) * ({r_s_frac[0]}/{r_s_frac[1]}) * ({rho_s_frac[0]}/{rho_s_frac[1]}) * 10^{total_exponent}")
    
    # 3. Simplify the combined fraction by canceling common factors
    mantissa_num, mantissa_den, simplified_nums, simplified_dens = simplify_fractions(numerators, denominators)

    print("\nSimplifying by canceling common factors across all fractions:")
    simplified_expr = " * ".join([f"({n}/{d})" for n, d in zip(simplified_nums, simplified_dens)])
    print(f"Simplified mantissa expression: {simplified_expr}")
    
    print(f"Final mantissa = {mantissa_num} / {mantissa_den} = {mantissa_num / mantissa_den}")
    
    # 4. Calculate the final force
    titan_force = (mantissa_num / mantissa_den) * (10**total_exponent)
    print(f"\nFinal Titan Force Calculation: F = {mantissa_num / mantissa_den} * 10^{total_exponent} = {titan_force} N")

    # 5. Calculate the true force using the full original formula for error analysis
    G_true = 6.6743e-11
    pi_true = math.pi
    m_p_true = 30.0
    r_c_true = 50e3  # 50 km in meters
    r_s_true = 1000e3 # 1000 km in meters
    rho_c_true = 1200.0 # 1.2 tons/m^3 in kg/m^3
    rho_s_true = 300.0  # 0.3 tons/m^3 in kg/m^3
    h_true = 500.0 # 500 m
    
    # Full mass calculation
    mass_core = (4/3) * pi_true * (r_c_true**3) * rho_c_true
    vol_shell = (4/3) * pi_true * (r_s_true**3 - r_c_true**3)
    mass_shell = vol_shell * rho_s_true
    M_true = mass_core + mass_shell

    # Full distance calculation
    r_true = r_s_true + h_true
    
    F_true = (G_true * M_true * m_p_true) / (r_true**2)
    
    # 6. Calculate absolute error and print the final answer
    abs_error = abs(F_true - titan_force)
    
    print(f"\nFor comparison, the true force is: {F_true:.4f} N")
    print(f"The absolute error is |{F_true:.4f} - {titan_force}| = {abs_error:.4f} N")
    
    # Final answer format
    print(f"\n<<<\Y{abs_error:.3f}>>>")

if __name__ == '__main__':
    main()