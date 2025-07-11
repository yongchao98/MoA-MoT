import math

def calculate_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Step 1: Define Given Information from the Plate Bearing Tests ---
    P_test_kN = 30.0  # Load during test in kN
    d_test_mm = 305.0 # Plate diameter in mm
    h_test_mm = 300.0 # Trial section thickness in mm
    delta_subgrade_um = 2460.0 # Deflection on subgrade in micrometers
    delta_pavement_um = 1080.0 # Deflection on trial section in micrometers

    # Convert units
    delta_subgrade_mm = delta_subgrade_um / 1000.0
    delta_pavement_mm = delta_pavement_um / 1000.0
    a_test_mm = d_test_mm / 2.0

    print("--- Plate Bearing Test Analysis ---")
    print(f"Test Plate Radius (a_test): {a_test_mm:.1f} mm")
    print(f"Trial Pavement Thickness (h_test): {h_test_mm:.1f} mm")
    print(f"Test Geometry Ratio (h_test / a_test): {h_test_mm / a_test_mm:.3f}")

    # For the same test setup, the influence factor F is the ratio of deflections.
    F_test = delta_pavement_mm / delta_subgrade_mm
    print(f"Settlement Influence Factor from Test (F_test): {F_test:.4f}")
    print("-" * 35)

    # --- Step 2: Define and Calculate Design Load Parameters ---
    W_design_ton = 1.80  # Design wheel load in metric tons
    p_design_kPa = 600.0 # Design tyre pressure in kN/m^2
    delta_allowable_mm = 1.00 # Allowable surface deflection in mm
    g = 9.81 # Acceleration due to gravity in m/s^2

    # Convert units
    P_design_N = W_design_ton * 1000 * g # Convert tons to kg, then to Newtons
    p_design_MPa = p_design_kPa / 1000.0 # Convert kN/m^2 to N/mm^2 (MPa)

    # Calculate radius of the design load area
    a_design_mm = math.sqrt(P_design_N / (math.pi * p_design_MPa))

    print("--- Design Load Analysis ---")
    print(f"Design Load (P_design): {P_design_N:.2f} N")
    print(f"Design Pressure (p_design): {p_design_MPa:.2f} N/mm^2")
    print(f"Design Load Radius (a_design): {a_design_mm:.2f} mm")
    print("-" * 35)
    
    # --- Step 3: Determine Required Pavement Thickness ---
    # The core principle is that for the same material properties (E1/E2),
    # the settlement factor F is a function of the h/a ratio.
    # If we assume the required F_design is the same as F_test, then the
    # h/a ratio must also be the same.
    # (h_design / a_design) = (h_test / a_test)
    # This allows us to solve for h_design directly.

    h_design_mm = h_test_mm * (a_design_mm / a_test_mm)

    print("--- Final Calculation ---")
    print("The required pavement thickness (h) is found by keeping the h/a ratio constant:")
    print("h_design / a_design = h_test / a_test")
    print("h_design = h_test * (a_design / a_test)\n")
    print("Substituting the values:")
    # Final output showing the equation with numbers
    print(f"h_design = {h_test_mm:.0f} mm * ({a_design_mm:.2f} mm / {a_test_mm:.1f} mm) = {h_design_mm:.2f} mm")


# Run the calculation
calculate_pavement_thickness()
# The final result rounded to one decimal place as is common in engineering.
final_answer = 300.0 * (math.sqrt(1.80 * 1000 * 9.81 / (math.pi * 0.6)) / (305.0 / 2.0))
print(f"<<<{final_answer:.1f}>>>")