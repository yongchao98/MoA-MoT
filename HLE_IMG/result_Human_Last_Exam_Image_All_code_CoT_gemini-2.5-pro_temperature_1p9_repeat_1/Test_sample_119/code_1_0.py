def determine_iupac_name():
    """
    This script outlines the step-by-step deduction of the compound's structure
    and IUPAC name based on the provided spectroscopic data.
    """

    # --- Step 1: Molecular Formula and Degree of Unsaturation ---
    print("--- Analysis Summary ---")
    print("\n1. Molecular Formula & DBE:")
    formula = "C9H13N"
    molecular_weight = 135
    C, H, N = 9, 13, 1
    
    print(f"   - Mass Spec M+ peak is {molecular_weight}, suggesting formula {formula}.")
    print(f"   - MW Check: {C}*12 + {H}*1 + {N}*14 = {C*12 + H*1 + N*14}")
    
    dbe = C - H/2 + N/2 + 1
    print(f"   - Degree of Unsaturation (DBE) = {C} - {H}/2 + {N}/2 + 1 = {dbe}. A DBE of 4 indicates a benzene ring.")

    # --- Step 2: Fragment Identification from NMR ---
    print("\n2. Fragment Identification:")
    c13_shifts = {'aromatic': [145.1, 128.5, 127.3, 126.3], 'aliphatic': [49.6, 43.5, 19.2]}
    h1_fragments = {
        'C6H5-': 'Multiplet, 5H, ~7.2 ppm',
        '-CH3': 'Doublet, 3H, ~1.2 ppm',
        '-CH-': 'Multiplet, 1H, ~2.8 ppm',
        '-CH2-': 'Multiplet, 2H, ~2.7 ppm',
        '-NH2': '2H (to satisfy formula C9H13N)'
    }
    
    print("   - From 13C NMR: A monosubstituted benzene ring and a 3-carbon side chain.")
    print("   - From DEPT-135 and 1H NMR, the side chain contains: one -CH3, one -CH2-, and one -CH- group.")
    print("   - The molecular formula requires an -NH2 group.")

    # --- Step 3: Structure Assembly ---
    print("\n3. Structure Assembly:")
    print("   - The 3H doublet at 1.2 ppm must be a -CH3 group next to a -CH group.")
    print("   - This gives the core fragment: -CH(CH3)-")
    print("   - The remaining fragments (C6H5-, -CH2-, -NH2) must be attached.")
    print("   - HSQC correlates the C/H signals:")
    print("     - C at 49.6 ppm is the -CH- group.")
    print("     - C at 43.5 ppm is the -CH2- group.")
    print("   - A C-N bond (~50 ppm) is more deshielded than a C-Ph bond (~40-45 ppm).")
    print("   - Therefore, the -NH2 is on the CH and the C6H5- is on the CH2.")
    print("   - Assembled Structure: C6H5-CH2-CH(NH2)-CH3")
    
    # --- Step 4: Final IUPAC Name ---
    iupac_name = "1-phenylpropan-2-amine"
    print("\n--- FINAL CONCLUSION ---")
    print(f"The IUPAC name of the compound is: {iupac_name}")

determine_iupac_name()