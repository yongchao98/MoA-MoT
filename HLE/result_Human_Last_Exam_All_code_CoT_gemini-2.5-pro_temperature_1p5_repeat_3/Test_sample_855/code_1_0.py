import math

def get_ground_truth_force():
    """Calculates the force using standard floating-point arithmetic for a baseline."""
    # Constants
    G = 6.674e-11  # Gravitational constant
    m_p = 50.0  # Probe mass in kg
    v_i = 300.0  # Initial velocity in m/s
    h = 5000.0  # Landing activation height in m

    # Pandora's properties
    rho_c = 1200.0  # Core density in kg/m^3
    rho_s = 300.0   # Shell density in kg/m^3
    r_c = 100e3     # Core radius in m
    a = 2000e3      # Equatorial radius in m
    c = 1985e3      # Polar radius in m

    # 1. Calculate net acceleration for landing
    # vf^2 = vi^2 + 2*a*d  => 0 = vi^2 + 2*a*h => a = -vi^2 / (2h)
    a_net = v_i**2 / (2 * h)  # Magnitude of upward acceleration

    # 2. Calculate Pandora's mass
    # V_spheroid = 4/3 * pi * a^2 * c
    # V_core = 4/3 * pi * r_c^3
    V_total = (4/3) * math.pi * (a**2) * c
    V_core = (4/3) * math.pi * r_c**3
    
    # M_total = M_core + M_shell
    M_core = rho_c * V_core
    M_shell = rho_s * (V_total - V_core)
    M_total = M_core + M_shell
    
    # 3. Calculate gravitational acceleration 'g' at the equator (landing site)
    # g = G * M_total / a^2
    g_pandora = G * M_total / a**2
    
    # 4. Calculate total required force
    # F_net = F_rocket - F_gravity => m*a_net = F_rocket - m*g
    # F_rocket = m * (a_net + g)
    F_decel = m_p * a_net
    F_gravity = m_p * g_pandora
    F_total = F_decel + F_gravity
    
    return F_total, F_decel, F_gravity, g_pandora

def find_best_titan_approximation(target_value):
    """
    Finds the best fraction n/d and exponent e to represent target_value.
    Representation: (n/d) * 10^e, where n, d are in [0, 31].
    """
    min_abs_error = float('inf')
    best_representation = (0, 1, 0) # n, d, e

    # We search for the best mantissa for exponents that are plausible for our result
    # The force is ~458, so an exponent of 2 (for 10^2) is the most likely candidate.
    exponent = 2
    target_mantissa = target_value / (10**exponent)
    
    best_n, best_d = 0, 1
    min_mantissa_error = float('inf')

    # Search for the best fraction (n/d) to approximate the target mantissa
    for d in range(1, 32):
        for n in range(0, 32):
            # Titan registers are 0-31
            current_mantissa = n / d
            error = abs(current_mantissa - target_mantissa)
            if error < min_mantissa_error:
                min_mantissa_error = error
                best_n, best_d = n, d
    
    best_titan_value = (best_n / best_d) * (10**exponent)
    min_abs_error = abs(best_titan_value - target_value)
    
    return (best_n, best_d, exponent), best_titan_value, min_abs_error

def main():
    """
    Simulates the calculation on the Titan architecture and finds the smallest error.
    """
    print("[Titan Feasibility Analysis]")
    
    # Step 1: Calculate ground truth values
    ground_truth_force, f_decel, f_gravity, g_pandora = get_ground_truth_force()
    
    print(f"\n1. Precise Off-board Calculation:")
    print(f"   - Gravitational Acceleration (g): {g_pandora:.4f} m/s^2")
    print(f"   - Deceleration Force (m * a_net): {f_decel:.4f} N")
    print(f"   - Gravitational Force (m * g): {f_gravity:.4f} N")
    print(f"   - Total Required Force (Ground Truth): {ground_truth_force:.4f} N")

    print("\n2. Titan On-board Calculation Strategy:")
    
    # Step 2: Calculate deceleration force on Titan
    # F_decel = 50 * (300^2 / (2*5000)) = 50 * 9 = 450.
    # To represent 450 in Titan sci-notation: 450 = 4.5 * 10^2.
    # The mantissa 4.5 can be represented exactly as the fraction 9/2.
    # Numerator 9 and denominator 2 are both valid 5-bit integers.
    f_decel_titan = "(9/2) * 10^2"
    print(f"   - Deceleration force is 450 N. This can be represented as {f_decel_titan}.")

    # Step 3: Calculate gravitational force on Titan
    # Calculating g from first principles on Titan would lead to massive overflows.
    # Instead, we approximate g with a simple, accurate fraction.
    # g = 0.1657 m/s^2. 1/6 = 0.1666... This is an excellent 5-bit fraction approximation.
    g_frac_n, g_frac_d = 1, 6
    print(f"   - Gravitational acceleration g ({g_pandora:.4f}) is best approximated by the fraction {g_frac_n}/{g_frac_d}.")
    # F_gravity = m * g = 50 * (1/6) = 50/6 = 25/3.
    # Numerator 25 and denominator 3 are valid 5-bit integers.
    f_gravity_titan_n, f_gravity_titan_d = 25, 3
    f_gravity_titan_val = f_gravity_titan_n / f_gravity_titan_d
    print(f"   - Gravitational force is {f_gravity_titan_val:.4f} N. Represented as the fraction {f_gravity_titan_n}/{f_gravity_titan_d}.")
    
    # Step 4: Final Summation and Representation
    # The final force value must be approximated as a single Titan representation.
    # F_total = 450 + 8.333... = 458.333... N
    # We find the best Titan scientific notation number to represent this value.
    print("\n3. Final Result Approximation:")
    print(f"   - The precise sum of forces is {ground_truth_force:.4f} N.")
    
    best_repr, best_val, min_error = find_best_titan_approximation(ground_truth_force)
    n, d, e = best_repr
    
    print(f"   - The best Titan approximation for this value is ({n}/{d}) * 10^{e}, which equals {best_val:.4f} N.")
    
    print("\n[Final Equation and Error]")
    # The problem asks to output the numbers in the final equation.
    print("The final calculated force is represented by the equation:")
    print(f"{n} / {d} * 10^{e} = {best_val:.1f} N")

    print(f"\nThe smallest absolute error achievable is {min_error:.3f} N.")
    
    # Output the final answer in the specified format
    final_answer = f"Y{min_error:.3f}"
    print(f"\nFinal Answer: {final_answer}")

if __name__ == "__main__":
    main()