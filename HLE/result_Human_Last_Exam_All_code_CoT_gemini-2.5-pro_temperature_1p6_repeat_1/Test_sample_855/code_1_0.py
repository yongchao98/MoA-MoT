import math

def display_titan_calculation():
    """
    This function demonstrates the calculation of the required rocket force
    using the constraints of the Titan 5-bit computer architecture.
    """

    # Step 1: Represent all input values in Titan's fractional scientific notation.
    # The format is (N/D) * 10^E, where N, D <= 31.
    m_probe_kg = "(1/2) * 10**2"  # 50 kg = 0.5 * 10^2
    v_initial_ms = "(3/1) * 10**2"  # 300 m/s = 3 * 10^2
    h_m = "(5/1) * 10**3"       # 5000 m = 5 * 10^3

    # Step 2: The required force is F_rocket = F_deceleration + F_gravity
    # F_rocket = m * a_decel + m * g
    
    # Step 3: Calculate the required deceleration, a_decel = v_i^2 / (2*h).
    # This calculation is possible within Titan's constraints.
    # v_i^2 = ((3/1) * 10**2)^2 = (9/1) * 10**4
    # 2*h = (2/1) * ((5/1) * 10**3) = (10/1) * 10**3 = (1/1) * 10**4
    # a_decel = ((9/1) * 10**4) / ((1/1) * 10**4) = (9/1) m/s^2
    a_decel_ms2 = "(9/1)"

    # Step 4: Calculate the force from deceleration, F_decel = m * a_decel.
    # F_decel = ((1/2) * 10**2) * (9/1) = (9/2) * 10**2 N = 450 N
    F_decel_N = "(9/2) * 10**2"

    # Step 5: Analyze Pandora's gravity, g.
    # A precise calculation of g ≈ 0.166 m/s^2 requires multiplying several fractional
    # approximations (for G, pi, r, rho), which results in intermediate numerators
    # (e.g., 20*21 = 420) that exceed the 5-bit limit of 31.
    #
    # We must therefore approximate or neglect g.
    # F_gravity = m*g = 50 * 0.166 ≈ 8.3 N.
    # F_deceleration = 450 N.
    # F_gravity is only about 1.8% of F_deceleration. Following the rule "eliminating
    # negligible terms," we can approximate g as 0 for the most accurate representable result.
    g_approx_ms2 = "0"
    
    # Step 6: The final force is approximated by neglecting gravity.
    # F_rocket ≈ F_decel
    
    print("--- Titan Landing Calculation ---")
    print("\n[Step 1: Input Values in Titan Format]")
    print(f"Probe Mass (m): {m_probe_kg} kg")
    print(f"Initial Velocity (v_i): {v_initial_ms} m/s")
    print(f"Altitude (h): {h_m} m")

    print("\n[Step 2: Kinematic Calculation]")
    print(f"Required deceleration (a) = v_i^2 / (2*h)")
    print(f"a = ({v_initial_ms})^2 / (2 * {h_m}) = {a_decel_ms2} m/s^2")

    print("\n[Step 3: Force Calculation & Simplification]")
    print("The gravitational force F_g (approx 8.3 N) is negligible compared to the deceleration force F_decel (450 N).")
    print("Per Titan rules, we eliminate the negligible term F_g (error ~1.8%).")
    
    print("\n[Step 4: Final Equation]")
    print(f"F_rocket ≈ m * a")
    print(f"F_rocket = ({m_probe_kg}) * ({a_decel_ms2})")
    
    final_result_str = f"F_rocket = {F_decel_N} N"
    print(final_result_str)
    
    # Step 7: Calculate the error of this approximation.
    # True value F_true ≈ 450 + 8.32 = 458.32 N.
    # Approximation F_approx = 450 N.
    # Absolute error e = |F_approx - F_true| / F_true
    e = abs(450 - 458.32) / 458.32
    # e_rounded = round(e, 3) # e = 0.01815 -> 0.018
    # The question requests the answer in a specific format outside the code block.

display_titan_calculation()

# The final answer format required by the user prompt.
# The smallest absolute error is |450 - 458.32| / 458.32 = 0.01815..., which rounds to 0.018.
print("\n<<<Y0.018>>>")