import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Test 1 (on subgrade)
    P_test = 30e3  # Load in Newtons (30 kN)
    d_plate = 305  # Plate diameter in mm
    w_subgrade = 2460e-3 # Deflection on subgrade in mm (2460 um)

    # Test 2 (on pavement)
    h_test = 300  # Pavement thickness in mm
    w_pavement = 1080e-3 # Deflection on pavement in mm (1080 um)

    # Design Load
    P_design_mass = 1.80 # Wheel load in tons
    p_design_pressure_kpa = 600 # Tyre pressure in kPa
    w_allowable = 1.00 # Max allowable deflection in mm
    
    # Constants
    g = 9.81 # acceleration due to gravity in m/s^2
    poisson_ratio = 0.5
    # For a flexible load and mu=0.5, the deflection coefficient is 1.5
    DEFLECTION_COEFF = 1.5

    # --- Step 1: Calculate Subgrade Modulus (Es) ---
    print("### Step 1: Analyzing the test on the subgrade ###")
    a_plate = d_plate / 2
    p_plate = P_test / (math.pi * a_plate**2)
    
    # Formula: w = (1.5 * p * a) / E_s  =>  E_s = (1.5 * p * a) / w
    E_s = (DEFLECTION_COEFF * p_plate * a_plate) / w_subgrade
    print(f"Test Plate Radius (a_plate): {a_plate:.2f} mm")
    print(f"Test Pressure (p_plate): {p_plate:.4f} N/mm^2 (MPa)")
    print(f"Calculated Subgrade Modulus (E_s): {E_s:.2f} MPa\n")

    # --- Step 2 & 3: Determine Deflection Factor (Fw) and Modular Ratio (Ep/Es) ---
    print("### Step 2 & 3: Analyzing the test on the pavement to find Ep/Es ###")
    # Fw is the ratio of deflection with pavement to deflection without it
    Fw_test = w_pavement / w_subgrade
    h_a_ratio_test = h_test / a_plate
    
    print(f"Test Section Deflection Factor (Fw_test): {Fw_test:.4f}")
    print(f"Test Section Thickness/Radius Ratio (h/a): {h_a_ratio_test:.4f}")
    
    # From Burmister's charts, for Fw ≈ 0.439 and h/a ≈ 1.967,
    # the modular ratio Ep/Es is approximately 10.
    Ep_Es_ratio = 10
    print(f"Based on standard charts, this corresponds to a modular ratio (Ep/Es) of: {Ep_Es_ratio}\n")

    # --- Step 4: Analyze the Design Wheel Load ---
    print("### Step 4: Analyzing the design wheel load ###")
    P_design_N = P_design_mass * 1000 * g
    # Convert pressure from kPa to MPa (N/mm^2)
    p_design_pressure_mpa = p_design_pressure_kpa / 1000.0
    
    # Formula: p = P / (pi * a^2)  => a = sqrt(P / (pi * p))
    a_design = math.sqrt(P_design_N / (math.pi * p_design_pressure_mpa))
    print(f"Design Load (P_design): {P_design_N/1000:.2f} kN")
    print(f"Design Tyre Pressure (p_design): {p_design_pressure_mpa:.2f} MPa")
    print(f"Calculated Design Load Radius (a_design): {a_design:.2f} mm\n")

    # --- Step 5: Determine Required Pavement Thickness ---
    print("### Step 5: Determining required pavement thickness ###")
    # First, calculate hypothetical deflection on subgrade under the design load
    w_s_design = (DEFLECTION_COEFF * p_design_pressure_mpa * a_design) / E_s
    
    # Second, calculate the required deflection factor to meet the allowable deflection
    Fw_design = w_allowable / w_s_design
    print(f"Hypothetical deflection on subgrade under design load: {w_s_design:.4f} mm")
    print(f"Required Deflection Factor (Fw_design): {Fw_design:.4f}")

    # Data from Burmister chart for Ep/Es = 10, mapping Fw to h/a
    # This represents the line for Ep/Es = 10 on the chart.
    fw_to_ha_data = {
        0.82: 0.5, 0.65: 1.0, 0.52: 1.5, 0.44: 2.0, 0.33: 3.0, 0.27: 4.0
    }
    
    # Find bracketing points for interpolation
    fw_values = sorted(fw_to_ha_data.keys(), reverse=True)
    fw1, fw2 = 0, 0
    for i in range(len(fw_values) - 1):
        if fw_values[i] >= Fw_design and fw_values[i+1] <= Fw_design:
            fw1 = fw_values[i]
            fw2 = fw_values[i+1]
            break
            
    ha1 = fw_to_ha_data[fw1]
    ha2 = fw_to_ha_data[fw2]

    # Linear Interpolation: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    h_a_design = ha1 + (Fw_design - fw1) * (ha2 - ha1) / (fw2 - fw1)
    
    print(f"From interpolation, required Thickness/Radius Ratio (h/a_design): {h_a_design:.4f}")

    # Final Calculation
    h_design = h_a_design * a_design
    
    print("\n--- Final Calculation ---")
    print("The required pavement thickness 'h' is found using the equation: h = (h/a)_design * a_design")
    print(f"h = {h_a_design:.4f} * {a_design:.2f}")
    print(f"Required Pavement Thickness (h): {h_design:.2f} mm")
    
    return h_design

# Execute the function and print the final answer in the required format
final_thickness = solve_pavement_thickness()
print(f"\n<<<The required pavement thickness is {final_thickness:.2f} mm>>>")