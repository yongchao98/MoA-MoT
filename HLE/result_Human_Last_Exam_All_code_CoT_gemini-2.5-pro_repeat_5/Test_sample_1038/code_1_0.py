import math

def calculate_resistance():
    """
    Calculates and prints the effective resistance of a subthreshold transistor
    with and without body biasing to demonstrate the benefit of Option C.
    """

    # --- Constants and Parameters ---
    n = 1.5      # Subthreshold slope factor
    U_T = 0.026  # Thermal voltage in Volts (approx. at 300K)
    # Technology-dependent current factor, approximates (W/L) * mu * Cox * U_T^2
    I_0 = 2e-7   # Amperes
    Vds = 0.05   # A small drain-source voltage in Volts for resistance calculation

    # --- Transistor Threshold Voltages ---
    Vt_normal = 0.45  # Normal threshold voltage in Volts
    Vt_biased = 0.65  # Threshold voltage with reverse body bias (increased Vt)

    # --- Operating Point ---
    # A gate-source voltage chosen to be below the normal Vt.
    Vgs_op = 0.4     # Operating gate-source voltage in Volts

    print("Analysis of Body Biasing on Pseudo-Resistor Performance:")
    print("-" * 60)
    print("This script calculates the resistance of a transistor in the subthreshold region.")
    print(f"The subthreshold current is modeled by the equation:")
    print("Ids = I_0 * exp((Vgs - Vt) / (n * U_T)) * (1 - exp(-Vds / U_T))")
    print(f"\nShared Parameters: Vgs={Vgs_op}V, Vds={Vds}V, n={n}, I_0={I_0:.1e}A")
    print("-" * 60)

    # --- Calculation for Normal Vt (Case 1) ---
    term1_exp = (Vgs_op - Vt_normal) / (n * U_T)
    term1_lin = (1 - math.exp(-Vds / U_T))
    Ids1 = I_0 * math.exp(term1_exp) * term1_lin
    R1 = Vds / Ids1

    print("\nCase 1: Standard Transistor (No Body Bias)")
    print(f"Inputs: Vt = {Vt_normal} V")
    print(f"Calculation for Resistance R1:")
    print(f"Ids = {I_0:.1e} * exp(({Vgs_op} - {Vt_normal}) / ({n} * {U_T})) * (1 - exp(-{Vds} / {U_T}))")
    print(f"Ids = {Ids1:.3e} A")
    print(f"R1 = {Vds} V / {Ids1:.3e} A = {R1 / 1e6:.2f} MOhms")

    # --- Calculation for Biased Vt (Case 2 - Option C) ---
    term2_exp = (Vgs_op - Vt_biased) / (n * U_T)
    term2_lin = (1 - math.exp(-Vds / U_T))
    Ids2 = I_0 * math.exp(term2_exp) * term2_lin
    R2 = Vds / Ids2

    print("\nCase 2: With Reverse Body Bias (Option C)")
    print(f"Inputs: Vt = {Vt_biased} V")
    print(f"Calculation for Resistance R2:")
    print(f"Ids = {I_0:.1e} * exp(({Vgs_op} - {Vt_biased}) / ({n} * {U_T})) * (1 - exp(-{Vds} / {U_T}))")
    print(f"Ids = {Ids2:.3e} A")
    print(f"R2 = {Vds} V / {Ids2:.3e} A = {R2 / 1e9:.2f} GOhms")

    print("-" * 60)
    print("\nConclusion:")
    print(f"Increasing Vt from {Vt_normal}V to {Vt_biased}V increases the effective resistance")
    print(f"from {R1/1e6:.2f} MOhms to {R2/1e9:.2f} GOhms, a factor of over {R2/R1:.0f}x.")
    print("This demonstrates how body biasing provides a much higher resistance,")
    print("which also leads to lower intrinsic leakage and better stability.")

calculate_resistance()

<<<C>>>