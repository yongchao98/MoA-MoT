import math

def solve_color_problem():
    """
    Calculates the properties of a Bromophenol blue solution to determine its color.
    """

    # --- Step 1: Define Constants and Initial Values ---
    # NOTE: The problem states 13 mg/ml, which is an extremely high concentration that would
    # result in an opaque solution. We assume it's a typo and should be 13 ug/ml.
    conc_ug_per_ml = 13.0  # Assumed concentration in µg/mL
    mw = 669.96  # Molar mass of Bromophenol blue in g/mol
    pKa = 4.0  # pKa of Bromophenol blue
    ph = 4.21  # pH of the solution

    # Path lengths in cm
    b_thin_mm = 1.0  # in mm
    b_thick_cm = 10.0 # in cm
    b_thin_cm = b_thin_mm / 10.0

    # Molar absorptivity (ε) values in M⁻¹cm⁻¹ at key wavelengths (approximate values)
    # HIn is the yellow Acidic form, In- is the blue Basic form
    # Wavelength 436nm (where yellow form absorbs best)
    epsilon_HIn_436 = 20000
    epsilon_In_436 = 10000
    # Wavelength 590nm (where blue form absorbs best)
    epsilon_HIn_590 = 1000
    epsilon_In_590 = 65000

    print("--- Step 1: Initial Parameters ---")
    print(f"Assumed Concentration: {conc_ug_per_ml} µg/mL")
    print(f"Molar Mass: {mw} g/mol")
    print(f"pH: {ph}, pKa: {pKa}")
    print(f"Path Lengths: {b_thin_cm} cm (thin) and {b_thick_cm} cm (thick)\n")

    # --- Step 2: Calculate Total Molar Concentration (Molarity) ---
    conc_g_per_L = conc_ug_per_ml / 1000.0  # convert µg/mL to g/L
    c_total = conc_g_per_L / mw
    print("--- Step 2: Total Molar Concentration ---")
    print(f"Total concentration = ({conc_g_per_L:.4f} g/L) / ({mw} g/mol) = {c_total:.4e} M\n")

    # --- Step 3: Calculate Ratio of Basic to Acidic Form ---
    # Using Henderson-Hasselbalch: pH = pKa + log([In-]/[HIn])
    # => log([In-]/[HIn]) = pH - pKa
    log_ratio = ph - pKa
    ratio = 10**log_ratio
    print("--- Step 3: Ratio of Indicator Forms ([Blue]/[Yellow]) ---")
    print(f"[Blue]/[Yellow] = 10^(pH - pKa) = 10^({ph} - {pKa}) = 10^{log_ratio:.2f} = {ratio:.4f}\n")

    # --- Step 4: Calculate Individual Concentrations ---
    # c_total = [HIn] + [In-] = [HIn] + ratio * [HIn] = [HIn] * (1 + ratio)
    c_HIn = c_total / (1 + ratio)  # Concentration of acidic form (Yellow)
    c_In = ratio * c_HIn           # Concentration of basic form (Blue)
    print("--- Step 4: Individual Concentrations of Forms ---")
    print(f"Concentration of Yellow form [HIn] = {c_HIn:.4e} M")
    print(f"Concentration of Blue form [In⁻] = {c_In:.4e} M\n")

    # --- Step 5: Calculate Absorbance (A = εbc) at Both Wavelengths and Path Lengths ---
    # Absorbance per cm of path length
    abs_per_cm_436 = (epsilon_HIn_436 * c_HIn) + (epsilon_In_436 * c_In)
    abs_per_cm_590 = (epsilon_HIn_590 * c_HIn) + (epsilon_In_590 * c_In)
    
    # Calculate total absorbance for thin path (b = 0.1 cm)
    A_thin_436 = abs_per_cm_436 * b_thin_cm
    A_thin_590 = abs_per_cm_590 * b_thin_cm

    # Calculate total absorbance for thick path (b = 10 cm)
    A_thick_436 = abs_per_cm_436 * b_thick_cm
    A_thick_590 = abs_per_cm_590 * b_thick_cm
    
    print("--- Step 5: Absorbance Calculations ---")
    print(f"Since both yellow and blue forms are present, the solution color is green.")
    print("\nFor the THIN side (b = 1.0 mm = 0.1 cm):")
    print(f"Absorbance at 436nm (absorbs blue light) = (({epsilon_HIn_436} * {c_HIn:.2e}) + ({epsilon_In_436} * {c_In:.2e})) * {b_thin_cm} = {A_thin_436:.4f}")
    print(f"Absorbance at 590nm (absorbs yellow light) = (({epsilon_HIn_590} * {c_HIn:.2e}) + ({epsilon_In_590} * {c_In:.2e})) * {b_thin_cm} = {A_thin_590:.4f}")
    print("Interpretation: Low absorbance values result in a LIGHT GREEN color.")

    print("\nFor the THICK side (b = 10 cm):")
    print(f"Absorbance at 436nm (absorbs blue light) = (({epsilon_HIn_436} * {c_HIn:.2e}) + ({epsilon_In_436} * {c_In:.2e})) * {b_thick_cm} = {A_thick_436:.4f}")
    print(f"Absorbance at 590nm (absorbs yellow light) = (({epsilon_HIn_590} * {c_HIn:.2e}) + ({epsilon_In_590} * {c_In:.2e})) * {b_thick_cm} = {A_thick_590:.4f}")
    print("Interpretation: High absorbance values result in a saturated, dark GREEN color.\n")

    # --- Step 6: Final Conclusion ---
    print("--- Conclusion ---")
    print("The color is green due to the mix of yellow and blue forms of Bromophenol blue.")
    print("Through the thin side (1mm), the path length is short, absorbance is low, and the color is light.")
    print("Through the thick side (10cm), the path length is long, absorbance is high, and the color is dark/saturated.")
    print("Result -> Thin: light green, Thick: green.")

solve_color_problem()
<<<C>>>