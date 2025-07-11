import math

def calculate_escape_velocity_on_titan():
    """
    This function simulates the calculation of Pandora's escape velocity
    using the constraints of the Titan 4-bit architecture.
    """

    # Step 1: Define constants and approximations
    # We simplify the problem by ignoring the core's higher density.
    # The error from this simplification is small.
    # v_e^2 = (8 * pi * G * R^2 * rho_s) / 3
    # We group the constants to manage the 4-bit constraint.
    # (8 * pi * G) / 3 has a true value of ~5.592e-10.
    # We approximate this with 11/2 = 5.5, which is a good fit and calculable.
    
    # Term: (8 * pi * G) / 3
    # Approximation: 11/2 * 10^-10
    term1_num = 11
    term1_den = 2
    term1_exp = -10

    # Term: R^2 = (2e6)^2 = 4e12
    # Mantissa for R^2
    term2_num = 4
    term2_exp = 12

    # Term: rho_s = 300 = 3e2
    # Mantissa for rho_s
    term3_num = 3
    term3_exp = 2

    # Step 2: Simulate Titan's multiplication for v_e^2 = term1 * term2 * term3
    # All intermediate numerators/denominators must be <= 15.

    # First multiplication: (11/2 * 10^-10) * (4 * 10^12)
    # Mantissa: 11/2 * 4 = 22. This is > 15.
    # Per Titan rules, we expand: 22 = 14 + 8.
    # This expression is held in a register.
    interim_mantissa_expr = [14, 8]
    interim_exp = term1_exp + term2_exp # -10 + 12 = 2

    # Second multiplication: (14+8)e2 * (3e2)
    # Mantissa: (14 * 3) + (8 * 3)
    # 14 * 3 = 42 > 15. Expand: 42 = 14 + 14 + 14
    # 8 * 3 = 24 > 15. Expand: 24 = 8 + 8 + 8
    final_mantissa_expr = [14, 14, 14, 8, 8, 8]
    final_exp = interim_exp + term3_exp # 2 + 2 = 4
    
    # The value for v_e^2 is (14+14+14+8+8+8) * 10^4 = 66 * 10^4
    v_sq_mantissa = sum(final_mantissa_expr)
    v_sq_exp = final_exp
    v_sq_titan = v_sq_mantissa * (10**v_sq_exp)

    # Step 3: Approximate the square root using one step of Newton-Raphson
    # v_e = sqrt(v_sq_mantissa) * 10^(v_sq_exp/2)
    # We need to calculate sqrt(66).
    # Formula: y_next = 0.5 * (y_prev + N / y_prev)
    # Initial guess y_prev for sqrt(66) is 8 (since 8*8=64).
    y_prev_num = 8
    
    # N / y_prev = 66 / 8 = 33 / 4 = 8 + 1/4. This is calculable.
    N_div_y_num = 8
    N_div_y_den = 4 # Represents 1/4
    
    # y_prev + N/y_prev = 8 + (8 + 1/4) = 16 + 1/4.
    # 16 > 15. We must handle this during the next step.
    
    # 0.5 * (16 + 1/4) = (0.5 * 16) + (0.5 * 1/4) = 8 + 1/8.
    # This is a valid result on Titan.
    y_next_whole = 8
    y_next_frac_num = 1
    y_next_frac_den = 8
    
    # Final calculated mantissa for v_e
    v_e_mantissa_titan = y_next_whole + y_next_frac_num / y_next_frac_den
    v_e_exp_titan = v_sq_exp / 2
    
    v_e_titan = v_e_mantissa_titan * (10**v_e_exp_titan)

    # Step 4: Calculate the "true" value for comparison
    G_true = 6.67430e-11
    PI_true = math.pi
    R_total_true = 2e6
    r_core_true = 1e5
    rho_core_true = 1200
    rho_shell_true = 300
    
    vol_core = (4/3) * PI_true * (r_core_true**3)
    mass_core = vol_core * rho_core_true
    
    vol_shell = (4/3) * PI_true * (R_total_true**3 - r_core_true**3)
    mass_shell = vol_shell * rho_shell_true
    
    M_true = mass_core + mass_shell
    
    v_e_true = math.sqrt((2 * G_true * M_true) / R_total_true)
    
    # Step 5: Calculate the absolute error
    error = abs(v_e_titan - v_e_true)

    # Step 6: Print the results as requested
    print("Calculation of Pandora's Escape Velocity on Titan Architecture")
    print("-" * 60)
    print("Formula: v_e = sqrt(2*G*M/R)")
    print("Simplified Formula: v_e ≈ sqrt( (8*pi*G/3) * R^2 * rho_s )")
    print("\nApproximations used:")
    print(f"  (8*pi*G/3) ≈ {term1_num}/{term1_den} * 10^{term1_exp}")
    print(f"  R^2 mantissa ≈ {term2_num}")
    print(f"  rho_s mantissa ≈ {term3_num}")
    
    print("\nFinal Equation on Titan:")
    # This line prints each number used in the final equation.
    print(f"v_e = sqrt( ( {term1_num}/{term1_den} * 10^{term1_exp} ) * ( {term2_num} * 10^{term2_exp} ) * ( {term3_num} * 10^{term3_exp} ) )")
    print(f"    ≈ ( {y_next_whole} + {y_next_frac_num}/{y_next_frac_den} ) * 10^{int(v_e_exp_titan)}")
    
    print("\nResults:")
    print(f"  Titan Calculated v_e: {v_e_titan:.2f} m/s")
    print(f"  High-Precision v_e: {v_e_true:.2f} m/s")
    print(f"  Absolute Error: {error:.2f} m/s")
    
    # Final answer in the required format
    print("\n" + "<<" + f"Y[{error:.2f}]" + ">>")

calculate_escape_velocity_on_titan()