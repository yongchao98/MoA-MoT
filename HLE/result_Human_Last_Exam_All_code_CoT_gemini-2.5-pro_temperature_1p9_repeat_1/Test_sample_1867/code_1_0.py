import math

def calculate_transformation_ratio():
    """
    Calculates the impedance transformation ratio for a diode signal source.
    """
    # --- Input Parameters ---
    V1 = 0.78         # Start voltage of linear region in Volts
    V2 = 0.98         # End voltage of linear region in Volts
    delta_I = 0.445   # Change in current in the linear region in Amps
    R_load = 50.0     # Load resistance in Ohms
    margin = 0.20     # Startup margin (20%)

    print("--- Diode Impedance Matching Calculation ---\n")

    # Step 1: Calculate the change in voltage (ΔV) across the linear region.
    delta_V = V2 - V1
    print(f"Step 1: Calculate the change in voltage (ΔV)")
    print(f"ΔV = {V2} V - {V1} V = {delta_V:.2f} V\n")

    # Step 2: Calculate the nominal dynamic resistance of the diode (the source resistance).
    # This is constant in the linear region: Rs_nominal = ΔV / ΔI.
    Rs_nominal = delta_V / delta_I
    print(f"Step 2: Calculate the nominal source resistance (Rs_nominal)")
    print(f"Rs_nominal = ΔV / ΔI = {delta_V:.2f} V / {delta_I} A = {Rs_nominal:.4f} Ohms\n")

    # Step 3: Apply the 20% startup margin to find the design source resistance.
    Rs_design = Rs_nominal * (1 + margin)
    print(f"Step 3: Apply the {margin*100:.0f}% startup margin to get the design resistance (Rs_design)")
    print(f"Rs_design = Rs_nominal * (1 + {margin}) = {Rs_nominal:.4f} Ohms * {1+margin} = {Rs_design:.4f} Ohms\n")

    # Step 4: Calculate the impedance transformation ratio for optimum power transfer.
    # For maximum power transfer, the transformed load impedance must match the source impedance.
    # Transformation Ratio (N) = R_source_design / R_load
    transformation_ratio = Rs_design / R_load
    print(f"Step 4: Calculate the final impedance transformation ratio (N)")
    print(f"For optimum power transfer, N = Rs_design / R_load")
    print(f"The final equation: N = ({Rs_design:.4f} Ohms) / ({R_load} Ohms)")
    print(f"N = {transformation_ratio:.5f}\n")

    print("------------------------------------------")
    print("The final impedance transformation ratio is:")
    print(f"{transformation_ratio:.5f}")


calculate_transformation_ratio()

# Final numeric answer provided separately for verification
# The required transformation ratio is ( (V2-V1) / ΔI ) * (1 + margin) / R_load
# ( (0.98 - 0.78) / 0.445 ) * 1.2 / 50 = 0.0107865...
# Rounded to five decimal places, the answer is 0.01079