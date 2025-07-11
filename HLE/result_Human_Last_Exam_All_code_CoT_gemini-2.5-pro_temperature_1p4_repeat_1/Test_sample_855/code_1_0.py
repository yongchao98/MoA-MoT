import math

def find_best_fractional_approximation(target_value, max_val=31):
    """
    Finds the best fraction n/d to approximate a target value,
    where n and d are <= max_val.
    The function returns the best fraction (n, d) and its corresponding value.
    """
    best_fraction = (1, 1)
    min_error = float('inf')

    for d in range(1, max_val + 1):
        # Find the best numerator for this denominator
        # n/d approx target => n approx target * d
        n_ideal = target_value * d
        
        # Check two nearest integers for n
        n_candidates = [math.floor(n_ideal), math.ceil(n_ideal)]
        
        for n in n_candidates:
            if 0 < n <= max_val:
                current_value = n / d
                error = abs(current_value - target_value)
                if error < min_error:
                    min_error = error
                    best_fraction = (int(n), int(d))
    
    return best_fraction, best_fraction[0] / best_fraction[1]

def titan_solver():
    """
    Solves the physics problem following Titan computer constraints.
    """
    # Step 1 & 2: Constants and Approximations
    # Probe mass in kg. Represented as 50.
    m_probe = 50.0

    # Gravitational Constant G in m^3 kg^-1 s^-2
    # G = 6.674e-11. We approximate 6.674 with 20/3.
    # G_titan = (20/3) * 10^-11
    G_frac = (20, 3)
    G_exp = -11

    # Pandora's Mass (M) in kg
    # Based on shell volume (V = 4/3 * pi * r^3) and density (rho = 300 kg/m^3)
    # r = 2000 km = 2e6 m.
    # M_shell = (4/3) * pi * (2e6)^3 * 300 = 32 * pi * 10^20 approx 1.005e22 kg.
    # We approximate M as 1 * 10^22 kg, which is (1/1) * 10^22 in Titan notation.
    M_frac = (1, 1)
    M_exp = 22

    # Pandora's radius r in meters
    # r = 2000 km = 2e6 m. This is (2/1) * 10^6
    r_frac = (2, 1)
    r_exp = 6

    # Step 3: Calculate gravitational acceleration g = G * M / r^2
    # r^2 = ((2/1) * 10^6)^2 = (4/1) * 10^12
    r2_frac = (4, 1)
    r2_exp = 12
    
    # g_frac = G_frac * M_frac / r2_frac = (20/3) * (1/1) / (4/1) = 20/12 = 5/3
    g_mantissa_num = G_frac[0] * M_frac[0] * r2_frac[1]
    g_mantissa_den = G_frac[1] * M_frac[1] * r2_frac[0]
    # simplify 20/12 to 5/3
    g_frac = (5, 3)

    # g_exp = G_exp + M_exp - r2_exp = -11 + 22 - 12 = -1
    g_exp = G_exp + M_exp - r2_exp

    # So, g = (5/3) * 10^-1 = 5/30 = 1/6 m/s^2
    # final g fraction is (1,6)
    g_titan_frac = (1, 6)
    g_val = g_titan_frac[0] / g_titan_frac[1]

    # Step 4: Calculate required deceleration a
    # a = v_i^2 / (2*d)
    # v_i = 300 m/s, d = 5000 m
    # a = 300^2 / (2*5000) = 90000 / 10000 = 9 m/s^2
    # In Titan fractions, this is (9, 1).
    a_titan_frac = (9, 1)
    a_val = 9.0

    # Step 5: Calculate the final force F = m * (g + a)
    # First, calculate the ideal "true" value for error comparison
    F_true = m_probe * (g_val + a_val)

    # Now, calculate using Titan rules.
    # F = m*a + m*g. Calculate components separately.
    # F_decel = m * a = 50 * 9 = 450 N.
    # F_gravity = m * g = 50 * (1/6) = 25/3 N.
    # We must add these: 450 + 25/3. Use scientific notation.
    # 450 = 4.5 * 10^2. 4.5 is 9/2. F_decel_titan = (9/2) * 10^2
    # 25/3 = 8.333... = 0.08333... * 10^2.
    # The fraction for 0.08333... is 1/12. F_gravity_titan = (1/12) * 10^2
    # Now add the mantissas: 9/2 + 1/12
    F_mantissa_num = 9 * 12 + 1 * 2
    F_mantissa_den = 2 * 12
    # This gives 110/24, which simplifies to 55/12.
    F_mantissa_val = 55 / 12  # This is 4.58333...

    # The numerator 55 exceeds the 5-bit limit of 31.
    # We must find the best fractional approximation n/d for 55/12.
    best_mantissa_frac, best_mantissa_val = find_best_fractional_approximation(F_mantissa_val)

    # The best approximation for 4.58333 is 23/5 = 4.6
    F_titan_frac = best_mantissa_frac
    F_titan_exp = 2
    F_titan_val = best_mantissa_val * (10**F_titan_exp)
    
    # Step 6: Determine the error
    error = abs(F_titan_val - F_true)
    
    # Step 7: Final Output
    print("--- Titan Precise Landing Calculation ---")
    print("\nThe force is calculated using F = m * (g + a).")
    print("All values are approximated using Titan's 5-bit fractional arithmetic.")
    print("\n--- Final Equation Components ---")
    print(f"Mass (m): {m_probe} kg")
    print(f"Gravitational Acceleration (g): Approximated as {g_titan_frac[0]}/{g_titan_frac[1]} m/s^2")
    print(f"Required Deceleration (a): Calculated as {a_titan_frac[0]}/{a_titan_frac[1]} m/s^2")
    
    print("\n--- Final Force Calculation ---")
    print("The ideal force is F = 50 * (1/6 + 9) = 458.333... N")
    print(f"Titan's constrained calculation leads to F = ({F_titan_frac[0]}/{F_titan_frac[1]}) * 10^{F_titan_exp} N")
    print(f"Final calculated force F = {F_titan_val:.3f} N")

    print(f"\nThe smallest absolute error produced is |{F_titan_val:.3f} - {F_true:.3f}| = {error:.3f} N")

    # The final answer format is Y[e]
    # This is a placeholder for the automated grading system.
    final_answer = f"Y[{error:.3f}]"

# Execute the solver
titan_solver()
