import math

def solve_color_problem():
    """
    Calculates the properties of a bromophenol blue solution to determine its color
    based on pH and viewing path length.
    """
    # --- Step 1: Define constants and initial values ---
    pH = 4.21
    pKa_bromophenol_blue = 4.1
    conc_mass_mg_ml = 13.0  # mg/mL
    molar_mass_g_mol = 669.96 # g/mol
    path_thin_mm = 1.0       # mm
    path_thick_cm = 10.0     # cm

    print("--- Initial Parameters ---")
    print(f"Solution pH: {pH}")
    print(f"Bromophenol Blue pKa: {pKa_bromophenol_blue}")
    print(f"Path Length (Thin): {path_thin_mm} mm")
    print(f"Path Length (Thick): {path_thick_cm} cm\n")

    # --- Step 2: Determine the color from pH ---
    # At a pH close to the pKa, both acidic (yellow) and basic (blue) forms of the indicator exist.
    # The mixture of yellow and blue results in a green color.
    print("--- Color Determination ---")
    if abs(pH - pKa_bromophenol_blue) < 1.0:
        print("The pH is close to the pKa, so both the yellow (acidic) and blue (basic) forms are present.")
        print("Resulting color: Yellow + Blue => Green\n")

    # --- Step 3: Use Henderson-Hasselbalch to find the ratio of forms ---
    # pH = pKa + log10([In-]/[HIn])
    ratio_In_over_HIn = 10**(pH - pKa_bromophenol_blue)

    # --- Step 4: Calculate Molar Concentration ---
    # Convert mass concentration to molar concentration
    conc_mass_g_l = conc_mass_mg_ml # 13 mg/mL is equivalent to 13 g/L
    total_molar_conc = conc_mass_g_l / molar_mass_g_mol

    # Calculate concentration of each form
    # total_conc = [HIn] + [In-] = [HIn] + ratio * [HIn] = [HIn] * (1 + ratio)
    conc_HIn = total_molar_conc / (1 + ratio_In_over_HIn) # Acidic (yellow) form
    conc_In = total_molar_conc - conc_HIn                # Basic (blue) form

    print("--- Concentration Calculations ---")
    print(f"Ratio of [Blue Form]/[Yellow Form]: {ratio_In_over_HIn:.2f}")
    print(f"Concentration of Yellow Form [HIn]: {conc_HIn:.4f} M")
    print(f"Concentration of Blue Form [In-]: {conc_In:.4f} M\n")


    # --- Step 5: Compare perceived color based on path length using Beer-Lambert Law (A = ebc) ---
    path_thin_cm = path_thin_mm / 10.0
    path_thick_cm = path_thick_cm

    print("--- Path Length & Intensity ---")
    print("According to the Beer-Lambert Law (A = ebc), absorbance (A) is directly proportional to path length (b).")
    print(f"Path length through the thick side ({path_thick_cm} cm) is {path_thick_cm / path_thin_cm:.0f} times longer than the thin side ({path_thin_cm} cm).")
    print("A longer path length results in higher absorbance and a more intensely perceived color.\n")

    # --- Step 6: Conclusion ---
    print("--- Conclusion ---")
    print("The solution will appear green when viewed through either side.")
    print("Through the thin side (1mm), the color will be 'light green' due to lower absorbance.")
    print("Through the thick side (10cm), the color will be a much more intense 'green' due to higher absorbance.")
    print("\nThis corresponds to the answer: Thin: light green, Thick: green")


solve_color_problem()
<<<C>>>