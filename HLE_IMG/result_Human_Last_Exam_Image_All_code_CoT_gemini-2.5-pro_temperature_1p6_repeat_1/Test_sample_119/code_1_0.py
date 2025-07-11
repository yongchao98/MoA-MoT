def identify_compound():
    """
    This script details the step-by-step analysis of the provided spectral
    data to determine the IUPAC name of the unknown compound.
    """
    
    # Data provided from the problem description
    molecular_ion_peak = 135
    c13_nmr_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    dept_135_signals = {"negative": 1, "positive": 5}
    
    print("--- Structure Elucidation Report ---")
    
    # Step 1: Mass Spectrum Analysis
    print("\n1. Mass Spectrum (MS) Analysis:")
    print(f"The molecular ion (M+) peak is observed at m/z = {molecular_ion_peak}.")
    print("Based on the Nitrogen Rule, an odd molecular weight suggests an odd number of nitrogen atoms.")
    # Based on other data (NMR suggesting a phenyl group C6H5, mass 77), we can deduce the formula.
    # 135 (Total) - 77 (C6H5) - 14 (N) = 44 (C3H8).
    # Proposed formula: C9H13N.
    C, H, N = 9, 13, 1
    dbe = C - H/2 + N/2 + 1
    print(f"The molecular formula is determined to be C{C}H{H}N{N}.")
    print(f"The Degree of Unsaturation is calculated as: {C} - ({H}/2) + ({N}/2) + 1 = {int(dbe)}.")
    print("A DBE of 4 corresponds to a benzene ring.")
    
    # Step 2: 13C and DEPT-135 NMR Analysis
    print("\n2. 13C and DEPT-135 NMR Analysis:")
    print(f"The 13C NMR shows 7 distinct carbon signals: {c13_nmr_shifts}.")
    print(f"Since the formula is C9H13N, this indicates molecular symmetry (specifically, in the phenyl ring).")
    print(f"DEPT-135 data indicates {dept_135_signals['negative']} CH2 group (negative signal), {dept_135_signals['positive']} CH/CH3 groups (positive signals), and 1 quaternary carbon (absent from DEPT).")
    print("Assignments based on chemical shifts and DEPT data:")
    print(f"- 145.1 ppm: Quaternary aromatic C (C_ipso).")
    print(f"- 128.5, 127.3, 126.3 ppm: Aromatic CH carbons.")
    print(f"- 49.6 ppm: Aliphatic CH group (positive DEPT signal).")
    print(f"- 43.5 ppm: Aliphatic CH2 group (negative DEPT signal).")
    print(f"- 19.2 ppm: Aliphatic CH3 group (positive DEPT signal, most upfield).")

    # Step 3: 1H NMR and HSQC Analysis
    print("\n3. 1H NMR and HSQC Analysis:")
    print("The 1H NMR spectrum shows three main regions for the non-labile protons:")
    print("- ~7.2 ppm: A multiplet integrating to 5H, characteristic of a monosubstituted phenyl group (C6H5-).")
    print("- ~2.8 ppm: A complex multiplet integrating to 3H.")
    print("- ~1.2 ppm: A doublet integrating to 3H.")
    print("This pattern suggests a C6H5-CH2-CH-CH3 structure.")
    print("The HSQC spectrum confirms the direct C-H connections:")
    print(f"- The 1H signal at ~1.2 ppm is connected to the 13C signal at {19.2} ppm (CH3).")
    print(f"- The 1H signals around ~2.8 ppm are connected to the 13C signals at {43.5} ppm (CH2) and {49.6} ppm (CH).")
    
    # Step 4: Final Structure and IUPAC Name
    print("\n4. Conclusion:")
    print("Combining all the fragments (C6H5, CH2, CH, CH3, and a required NH2 group to match the formula) leads to a single possible structure.")
    print("Structure: C6H5-CH2-CH(NH2)-CH3")
    print("\nThe IUPAC name for this compound is:")
    print("1-phenylpropan-2-amine")

identify_compound()