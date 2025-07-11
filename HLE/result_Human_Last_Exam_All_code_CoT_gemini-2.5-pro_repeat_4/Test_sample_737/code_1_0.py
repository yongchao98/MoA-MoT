import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Step 1: Define Inputs and Constants ---
    # Plate bearing test data
    P_plate_kN = 30.0
    d_plate_mm = 305.0
    delta_subgrade_um = 2460.0
    delta_trial_um = 1080.0
    h_trial_mm = 300.0

    # Design load data
    W_design_ton = 1.80
    g_ms2 = 9.81  # acceleration due to gravity
    p_design_kPa = 600.0
    delta_design_mm = 1.00

    # Material property
    mu = 0.5

    # Convert units for calculation
    P_plate_N = P_plate_kN * 1000
    a_plate_mm = d_plate_mm / 2
    delta_subgrade_mm = delta_subgrade_um / 1000
    delta_trial_mm = delta_trial_um / 1000
    p_design_MPa = p_design_kPa / 1000 # N/mm^2
    P_design_N = W_design_ton * 1000 * g_ms2

    print("--- Problem Setup ---")
    print(f"Subgrade test deflection (Δ_subgrade): {delta_subgrade_mm} mm")
    print(f"Trial pavement deflection (Δ_trial): {delta_trial_mm} mm")
    print(f"Design deflection limit (Δ_design): {delta_design_mm} mm")
    print("-" * 25)

    # --- Step 2: Calculate Subgrade Modulus (E₂) from Rigid Plate Test ---
    p_plate_MPa = P_plate_N / (math.pi * a_plate_mm**2)
    # Using formula for rigid plate: Δ = (1.18 * p * a) / E
    E2_MPa = (1.18 * p_plate_MPa * a_plate_mm) / delta_subgrade_mm
    print(f"Step 2: Calculated Subgrade Modulus (E₂): {E2_MPa:.2f} MPa")

    # --- Step 3: Determine Modulus Ratio (E₁/E₂) ---
    # F₂ is the ratio of deflections under identical loading
    F2_trial = delta_trial_mm / delta_subgrade_mm
    h_a_ratio_trial = h_trial_mm / a_plate_mm
    print(f"Step 3a: Trial section h/a ratio: {h_a_ratio_trial:.3f}")
    print(f"Step 3b: Trial section deflection factor (F₂_trial): {F2_trial:.3f}")

    # This step simulates looking up E₁/E₂ from a Burmister chart.
    # We use linear interpolation on known chart data points for mu=0.5.
    # Data points are for F₂ at h/a=2.0 (close to our trial h/a of 1.967)
    # E₁/E₂=10 -> F₂=0.53; E₁/E₂=20 -> F₂=0.43
    E_ratio_1, F2_1 = 10, 0.53
    E_ratio_2, F2_2 = 20, 0.43
    E1_E2_ratio = E_ratio_1 + (E_ratio_2 - E_ratio_1) * (F2_trial - F2_1) / (F2_2 - F2_1)
    print(f"Step 3c: Interpolated Modulus Ratio (E₁/E₂): {E1_E2_ratio:.2f}")

    # --- Step 4: Analyze Design Wheel Load ---
    # Area = Load / Pressure
    a_design_mm_sq = P_design_N / (math.pi * p_design_MPa)
    a_design_mm = math.sqrt(a_design_mm_sq)
    print(f"Step 4: Design wheel load radius (a_design): {a_design_mm:.1f} mm")

    # --- Step 5: Determine Required Deflection Factor (F₂_design) ---
    # Using formula for flexible load: Δ = (1.5 * p * a) / E₂ * F₂
    # Rearranging for F₂: F₂ = (Δ * E₂) / (1.5 * p * a)
    F2_design = (delta_design_mm * E2_MPa) / (1.5 * p_design_MPa * a_design_mm)
    print(f"Step 5: Required design deflection factor (F₂_design): {F2_design:.4f}")

    # --- Step 6: Calculate Required Pavement Thickness (h_design) ---
    # This step finds the h/a ratio for our E₁/E₂ and F₂_design.
    # It requires interpolating on the Burmister chart data again.
    
    # First, construct points on our specific E₁/E₂=19.1 curve by interpolating known data
    def interpolate(x1, y1, x2, y2, x):
        return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

    # Known chart data at h/a=2 and h/a=3
    F2_at_ha2_E10, F2_at_ha3_E10 = 0.53, 0.42
    F2_at_ha2_E20, F2_at_ha3_E20 = 0.43, 0.33

    # F₂ on our E₁/E₂ curve at h/a = 2
    F2_at_ha2 = interpolate(E_ratio_1, F2_at_ha2_E10, E_ratio_2, F2_at_ha2_E20, E1_E2_ratio)
    # F₂ on our E₁/E₂ curve at h/a = 3
    F2_at_ha3 = interpolate(E_ratio_1, F2_at_ha3_E10, E_ratio_2, F2_at_ha3_E20, E1_E2_ratio)

    # Now, interpolate along the h/a axis to find the h/a ratio for our F₂_design
    h_a_1, h_a_2 = 2.0, 3.0
    h_a_design_ratio = interpolate(F2_at_ha2, h_a_1, F2_at_ha3, h_a_2, F2_design)
    print(f"Step 6a: Required design h/a ratio: {h_a_design_ratio:.3f}")

    # Final calculation
    h_design_mm = h_a_design_ratio * a_design_mm
    print("-" * 25)
    print("--- Final Result ---")
    print(f"Final calculation: Required Thickness = ({h_a_design_ratio:.3f}) * {a_design_mm:.1f} mm")
    print(f"Required pavement thickness = {h_design_mm:.1f} mm")
    return h_design_mm

# Run the calculation and store the result
final_thickness = solve_pavement_thickness()
# Print the final answer in the specified format
print(f"\n<<< {final_thickness:.1f} >>>")
