import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """

    # --- Given Data ---
    # Plate Bearing Test 1 (Subgrade)
    delta_subgrade = 2460e-6  # m (2460 μm)
    P_plate = 30e3           # N (30 kN)
    d_plate = 0.305          # m (305 mm)

    # Plate Bearing Test 2 (Trial Pavement)
    delta_trial = 1080e-6    # m (1080 μm)
    h_trial = 0.300          # m (300 mm)

    # Design Load
    W_design_ton = 1.80      # ton (metric)
    p_design = 600e3         # N/m^2 (600 kN/m^2)
    delta_allowable = 1.00e-3 # m (1.00 mm)

    # Material Properties
    mu = 0.5                 # Poisson's Ratio

    # Constants
    g = 9.81                 # m/s^2, acceleration due to gravity
    pi = math.pi

    # Burmister's Two-Layer Chart Data (F-factor for mu=0.5)
    # Source: Standard Pavement Engineering Textbooks
    # Format: {E1/E2 ratio: {h/a ratio: F factor}}
    chart = {
        10: {0.5: 0.82, 1.0: 0.67, 2.0: 0.48, 3.0: 0.38, 4.0: 0.30},
        20: {0.5: 0.74, 1.0: 0.58, 2.0: 0.36, 3.0: 0.25, 4.0: 0.20}
    }
    
    # --- Helper function for linear interpolation ---
    def lerp(x, x1, y1, x2, y2):
        return y1 + (x - x1) * (y2 - y1) / (x2 - x1)

    print("### STEP 1: Determine Subgrade Modulus (E₂) ###")
    a_plate = d_plate / 2
    # Deflection formula for single layer (mu=0.5): delta = (1.5 * P) / (pi * a * E)
    # Rearranging for E₂:
    E2 = (1.5 * P_plate) / (pi * a_plate * delta_subgrade)
    print(f"Plate radius (a_plate) = {a_plate:.4f} m")
    print(f"Subgrade Modulus (E₂) = (1.5 * {P_plate:.0f} N) / (π * {a_plate:.4f} m * {delta_subgrade} m)")
    print(f"Subgrade Modulus (E₂) = {E2 / 1e6:.2f} MPa\n")

    print("### STEP 2: Determine Pavement to Subgrade Modulus Ratio (E₁/E₂) ###")
    # For a two-layer system: F = (delta * pi * a * E₂) / (1.5 * P)
    F_trial = (delta_trial * pi * a_plate * E2) / (1.5 * P_plate)
    h_a_ratio_trial = h_trial / a_plate
    print(f"Trial section h/a ratio = {h_trial:.3f} m / {a_plate:.4f} m = {h_a_ratio_trial:.3f}")
    print(f"Calculated Deflection Factor (F_trial) from test = {F_trial:.4f}")

    # Interpolate to find E₁/E₂
    # First, find F for h/a = h_a_ratio_trial on the E₁/E₂=10 and E₁/E₂=20 curves
    F_10 = lerp(h_a_ratio_trial, 1.0, chart[10][1.0], 2.0, chart[10][2.0])
    F_20 = lerp(h_a_ratio_trial, 1.0, chart[20][1.0], 2.0, chart[20][2.0])
    print(f"Interpolated F on E₁/E₂=10 curve at h/a={h_a_ratio_trial:.3f} is {F_10:.4f}")
    print(f"Interpolated F on E₁/E₂=20 curve at h/a={h_a_ratio_trial:.3f} is {F_20:.4f}")
    
    # Now, interpolate E₁/E₂ using F_trial
    modular_ratio = lerp(F_trial, F_10, 10, F_20, 20)
    E1 = modular_ratio * E2
    print(f"Interpolated Modulus Ratio (E₁/E₂) for F_trial={F_trial:.4f} is {modular_ratio:.2f}")
    print(f"Pavement Modulus (E₁) = {modular_ratio:.2f} * {E2 / 1e6:.2f} MPa = {E1 / 1e6:.2f} MPa\n")

    print("### STEP 3: Determine Required Pavement Thickness (h_design) ###")
    P_design = W_design_ton * 1000 * g
    # Area = P / p = pi * a^2  => a = sqrt(P / (p * pi))
    a_design = math.sqrt(P_design / (p_design * pi))
    print(f"Design Load (P_design) = {W_design_ton} ton * 1000 kg/ton * {g} m/s² = {P_design:.2f} N")
    print(f"Design Load Radius (a_design) = sqrt({P_design:.2f} N / ({p_design:.0f} Pa * π)) = {a_design:.4f} m")

    # Calculate the required F-factor for the design
    F_design = (delta_allowable * pi * a_design * E2) / (1.5 * P_design)
    print(f"Required Deflection Factor (F_design) for allowable deflection = {F_design:.4f}")

    # Interpolate to find the design h/a ratio
    # First, find h/a for F_design on the E₁/E₂=10 and E₁/E₂=20 curves (inverse interpolation)
    ha_10 = lerp(F_design, chart[10][2.0], 2.0, chart[10][3.0], 3.0)
    ha_20 = lerp(F_design, chart[20][1.0], 1.0, chart[20][2.0], 2.0)
    print(f"Interpolated h/a on E₁/E₂=10 curve for F={F_design:.4f} is {ha_10:.4f}")
    print(f"Interpolated h/a on E₁/E₂=20 curve for F={F_design:.4f} is {ha_20:.4f}")
    
    # Now, interpolate the final h/a ratio using the known modular_ratio
    h_a_design_ratio = lerp(modular_ratio, 10, ha_10, 20, ha_20)
    print(f"Interpolated Design h/a Ratio for E₁/E₂={modular_ratio:.2f} is {h_a_design_ratio:.4f}\n")

    print("### STEP 4: Final Calculation ###")
    h_design = h_a_design_ratio * a_design
    h_design_mm = h_design * 1000

    print("The required pavement thickness is calculated as:")
    print(f"h_design = (h/a)_design * a_design")
    print(f"h_design = {h_a_design_ratio:.4f} * {a_design:.4f} m = {h_design:.4f} m")
    print(f"\nFinal Answer: The required pavement thickness is {h_design_mm:.1f} mm.")
    
    return h_design_mm

# Execute the function and get the final answer
final_thickness = solve_pavement_thickness()
print(f"<<<{final_thickness:.1f}>>>")