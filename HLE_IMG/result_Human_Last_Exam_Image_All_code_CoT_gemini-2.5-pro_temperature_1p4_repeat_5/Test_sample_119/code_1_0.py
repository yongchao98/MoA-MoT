import textwrap

def analyze_compound_spectra():
    """
    Analyzes the provided spectral data to determine the IUPAC name of the compound.
    This function will print the step-by-step analysis.
    """

    # --- Step 1: Mass Spectrometry Analysis ---
    print("--- Step 1: Mass Spectrometry (MS) Analysis ---")
    mw = 135
    print(f"The mass spectrum shows a molecular ion peak (M+) at m/z = {mw}.")
    print("According to the Nitrogen Rule, an odd molecular weight suggests the presence of an odd number of nitrogen atoms.")
    print("Let's assume one nitrogen atom (N, atomic weight ~14).")
    remaining_mass = mw - 14
    print(f"Remaining mass for Carbon and Hydrogen = {mw} - 14 = {remaining_mass}.")
    print("From other spectra (NMR), we see an aromatic ring (at least 6 carbons) and an alkyl chain.")
    print("A plausible molecular formula that fits this mass and contains a phenyl group is C9H13N.")
    print("Let's check the mass of C9H13N: (9 * 12) + (13 * 1) + (1 * 14) = 108 + 13 + 14 = 135. This matches.")
    molecular_formula = "C9H13N"
    print(f"Proposed Molecular Formula: {molecular_formula}\n")

    # --- Step 2: Degree of Unsaturation (DBE) ---
    print("--- Step 2: Degree of Unsaturation (DBE) Calculation ---")
    C, H, N = 9, 13, 1
    dbe = C + 1 - (H / 2) + (N / 2)
    print("The formula for DBE is: C + 1 - (H / 2) + (N / 2)")
    print(f"DBE = {C} + 1 - ({H} / 2) + ({N} / 2) = {int(dbe)}")
    print("A DBE of 4 is characteristic of a substituted benzene ring.\n")

    # --- Step 3: Infrared (IR) Spectroscopy Analysis ---
    print("--- Step 3: Infrared (IR) Spectroscopy Analysis ---")
    print(textwrap.dedent("""
    - Two sharp peaks around 3400-3300 cm⁻¹: Characteristic of a primary amine (R-NH2) N-H stretch.
    - Peaks just above 3000 cm⁻¹: Aromatic C-H stretch.
    - Peaks just below 3000 cm⁻¹: Aliphatic C-H stretch.
    - Peaks around 1600 & 1450 cm⁻¹: Aromatic C=C ring stretch.
    Conclusion from IR: The molecule contains a primary amine (-NH2), an aromatic ring, and an aliphatic part.\n"""))

    # --- Step 4: 13C NMR and DEPT-135 Analysis ---
    print("--- Step 4: 13C NMR and DEPT-135 Analysis ---")
    c13_shifts = {'145.1': 'Quaternary (C)', '128.5': 'CH', '127.3': 'CH', '126.3': 'CH',
                    '49.6': 'CH', '43.5': 'CH2 (negative)', '19.2': 'CH3'}
    print("13C Signals (ppm) and DEPT-135 assignment:")
    for shift, type in c13_shifts.items():
        print(f"- {shift}: {type}")
    print("\nAnalysis:")
    print("- Aromatic region (126-146 ppm): 4 signals suggest a monosubstituted benzene ring (1 ipso-C, 3 types of CH).")
    print("- Aliphatic region (19-50 ppm): 1 CH3, 1 CH2, and 1 CH group confirm a 3-carbon chain.\n")

    # --- Step 5: 1H NMR and HSQC Analysis ---
    print("--- Step 5: 1H NMR and HSQC Analysis ---")
    print(textwrap.dedent("""
    - ~7.3 ppm (5H, multiplet): A classic pattern for a monosubstituted benzene ring (C6H5-). HSQC confirms these protons are attached to the aromatic CH carbons.
    - ~2.8-2.9 ppm (total 3H, overlapping multiplets): HSQC shows these protons are on two different carbons.
        - One multiplet (1H) correlates to the CH carbon at 49.6 ppm.
        - Another multiplet (2H) correlates to the CH2 carbon at 43.5 ppm.
    - ~1.2 ppm (3H, doublet): This signal correlates to the CH3 carbon at 19.2 ppm. The doublet splitting indicates it is next to a single proton (a CH group).\n"""))
    
    # --- Step 6: Structure Assembly and Naming ---
    print("--- Step 6: Structure Assembly and Naming ---")
    print("Combining the pieces:")
    print("1. We have a phenyl group (C6H5-).")
    print("2. We have a primary amine (-NH2).")
    print("3. We have a propyl chain backbone: -CH-, -CH2-, -CH3.")
    print("\nPutting it all together:")
    print("- The CH3 group (~1.2 ppm) is a doublet, so it must be next to the CH group: -CH(R)-CH3.")
    print("- This leaves the CH2 group. The most logical structure is C6H5-CH2-CH(NH2)-CH3.")
    print("  - The CH2 group is benzylic (next to phenyl), explaining its shift (protons at ~2.8 ppm, carbon at 43.5 ppm).")
    print("  - The CH group is attached to the amine, the CH2, and the CH3, explaining its complex splitting and shift (proton at ~2.85 ppm, carbon at 49.6 ppm).")
    print("\nThe final structure is consistent with all spectral data.")
    final_name = "1-phenylpropan-2-amine"
    print(f"\nFinal IUPAC Name: {final_name}")

if __name__ == '__main__':
    analyze_compound_spectra()