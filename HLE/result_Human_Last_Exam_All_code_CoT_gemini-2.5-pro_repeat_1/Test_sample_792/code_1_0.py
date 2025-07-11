import math

def titan_calculation():
    """
    This function demonstrates the step-by-step calculation of the force F
    using the Titan computer's 5-bit fractional arithmetic rules.
    It prints the derivation and the final answer.
    """
    
    # --- Step 1: Define constants and approximations as 5-bit fractions ---
    # Physical parameters
    h_num, h_den = 10, 1  # h = 10m
    d_num, d_den = 20, 1  # d = 20m
    r_val_cm = 0.5 # r = 0.5 cm
    r_num, r_den = 1, 2   # r = 1/2 cm
    rho_val_kg_cm3 = 0.9 # rho = 0.9 kg/cm^3
    rho_num, rho_den = 9, 10 # rho = 9/10 kg/cm^3
    
    # Approximations for physical constants to fit 5-bit constraints
    # We choose simple values to ensure intermediate calculations are valid.
    # pi is approximated as 3/1. The real value is ~3.14159.
    pi_num, pi_den = 3, 1
    # g is approximated as 10/1. The real value is ~9.8 m/s^2.
    g_num, g_den = 10, 1

    print("Titan Calculation for Force F")
    print("="*35)
    print("F = (m * g * d) / h")
    print("where m = (4/3) * pi * r^3 * rho")
    print("\nStep 1: Define all values as 5-bit fractions.")
    print(f"h = {h_num}/{h_den} m")
    print(f"d = {d_num}/{d_den} m")
    print(f"r = {r_num}/{r_den} cm")
    print(f"rho = {rho_num}/{rho_den} kg/cm^3")
    print(f"pi ≈ {pi_num}/{pi_den}")
    print(f"g ≈ {g_num}/{g_den} m/s^2")

    # --- Step 2: Calculate mass (m) in kg ---
    # m = (4/3) * pi * (r)^3 * rho
    # We perform the calculation step-by-step to respect Titan's constraints.
    print("\nStep 2: Calculate mass (m).")
    print(f"m = (4/3) * ({pi_num}/{pi_den}) * ({r_num}/{r_den})^3 * ({rho_num}/{rho_den})")
    
    # First term: (4/3) * pi
    # (4/3) * (3/1) = 12/3 = 4/1. (Valid: 4 <= 31, 1 <= 31)
    m_calc1_num, m_calc1_den = 4 * pi_num, 3 * pi_den
    m_calc1_num, m_calc1_den = 4, 1 # Simplified
    print(f"  (4/3) * ({pi_num}/{pi_den}) = {m_calc1_num}/{m_calc1_den}")
    
    # Second term: m_calc1 * r^3
    # r^3 = (1/2)^3 = 1/8.
    # (4/1) * (1/8) = 4/8 = 1/2. (Valid: 1 <= 31, 2 <= 31)
    r3_num, r3_den = r_num**3, r_den**3
    m_calc2_num, m_calc2_den = m_calc1_num * r3_num, m_calc1_den * r3_den
    m_calc2_num, m_calc2_den = 1, 2 # Simplified
    print(f"  ({m_calc1_num}/{m_calc1_den}) * ({r_num}/{r_den})^3 = ({m_calc2_num}/{m_calc2_den})")

    # Final term for mass: m_calc2 * rho
    # (1/2) * (9/10) = 9/20. (Valid: 9 <= 31, 20 <= 31)
    m_num, m_den = m_calc2_num * rho_num, m_calc2_den * rho_den
    print(f"  ({m_calc2_num}/{m_calc2_den}) * ({rho_num}/{rho_den}) = {m_num}/{m_den}")
    print(f"Result: m = {m_num}/{m_den} kg")

    # --- Step 3: Calculate Force (F) ---
    print("\nStep 3: Calculate Force (F).")
    # Rearrange for safer calculation: F = (m * d / h) * g
    print(f"F = (({m_num}/{m_den}) * ({d_num}/{d_den}) / ({h_num}/{h_den})) * ({g_num}/{g_den})")
    
    # First term: m * d
    # (9/20) * (20/1) = 180/20. (Invalid: 180 > 31). Must simplify first.
    # (9/20) * (20/1) = 9/1 * 20/20 = 9/1. (Valid)
    f_calc1_num, f_calc1_den = m_num, m_den * (d_den/d_num) # Simplified
    f_calc1_num, f_calc1_den = 9, 1
    print(f"  ({m_num}/{m_den}) * ({d_num}/{d_den}) = {f_calc1_num}/{f_calc1_den}")
    
    # Second term: (m*d) / h
    # (9/1) / (10/1) = 9/10. (Valid)
    f_calc2_num, f_calc2_den = f_calc1_num * h_den, f_calc1_den * h_num
    print(f"  ({f_calc1_num}/{f_calc1_den}) / ({h_num}/{h_den}) = {f_calc2_num}/{f_calc2_den}")

    # Final term for force: (m*d/h) * g
    # (9/10) * (10/1) = 90/10. (Invalid: 90 > 31). Must simplify first.
    # (9/10) * (10/1) = 9/1. (Valid)
    F_num, F_den = f_calc2_num * g_num, f_calc2_den * g_den
    F_num, F_den = 9, 1 # Simplified
    print(f"  ({f_calc2_num}/{f_calc2_den}) * ({g_num}/{g_den}) = {F_num}/{F_den}")
    
    print("\n" + "="*35)
    print(f"Final Titan Calculation: F = {F_num}/{F_den} N")
    print("="*35)

    # --- Step 4: Calculate "Real" Force and Error ---
    m_real = (4/3) * math.pi * (r_val_cm**3) * (rho_val_kg_cm3 * 1e6) # convert rho to kg/m^3
    F_real = (m_real * 9.8 * d_num) / h_num
    error = abs(F_real - (F_num/F_den))
    
    print("\nFor Error Calculation:")
    print(f"Real-world mass m ≈ {m_real:.4f} kg")
    print(f"Real-world force F ≈ {F_real:.4f} N")
    print(f"Absolute error e = |{F_real:.4f} - {F_num}| = {error:.4f} N")

titan_calculation()