def solve_structure():
    """
    This script performs a step-by-step analysis of the provided spectroscopic data
    to determine the IUPAC name of the unknown compound.
    """
    # --- Provided Data ---
    molecular_weight = 135
    c13_shifts = [145.1, 128.5, 127.3, 126.3, 49.6, 43.5, 19.2]
    dept135_negative = 1
    dept135_positive = 5
    h1_shifts = {
        'aromatic': {'ppm': '~7.25', 'integration': 5, 'multiplicity': 'm'},
        'signal_A': {'ppm': '~2.8', 'integration': 2, 'multiplicity': 'm'},
        'signal_B': {'ppm': '~2.4', 'integration': 1, 'multiplicity': 'm'},
        'signal_C': {'ppm': '~1.2', 'integration': 3, 'multiplicity': 'd'}
    }
    hsqc_correlations = {
        '1.2': 19.2,
        '2.4': 49.6,
        '2.8': 43.5,
    }

    print("--- Step-by-Step Structure Elucidation ---")

    # --- Step 1: Molecular Formula from Mass Spectrometry ---
    print("\nStep 1: Molecular Formula Analysis")
    print(f"The molecular ion peak (M+) is at m/z = {molecular_weight}.")
    print("An odd molecular weight suggests an odd number of nitrogen atoms (Nitrogen Rule).")
    molecular_formula = "C9H13N"
    num_C = 9
    num_H = 13
    num_N = 1
    print(f"The molecular formula is determined to be {molecular_formula}.")
    
    # --- Step 2: Degree of Unsaturation (DBE) Calculation ---
    print("\nStep 2: Degree of Unsaturation (DBE)")
    dbe = num_C + 1 - (num_H / 2) + (num_N / 2)
    print(f"The DBE is calculated using the formula: C + 1 - H/2 + N/2")
    print(f"DBE = {num_C} + 1 - ({num_H}/2) + ({num_N}/2) = {int(dbe)}")
    print("A DBE of 4 strongly suggests a benzene ring.")

    # --- Step 3: NMR and Fragment Analysis ---
    print("\nStep 3: NMR Analysis and Structure Assembly")
    print("1H NMR shows a monosubstituted phenyl group (5H at ~7.25 ppm) and an aliphatic chain.")
    print("The aliphatic signals are a 3H doublet, a 1H multiplet, and a 2H multiplet.")
    print("This corresponds to a -CH3, a -CH, and a -CH2 group, respectively.")
    print("The coupling pattern (-CH3 is a doublet) indicates a -CH-CH3 fragment.")
    print("The full aliphatic backbone is therefore -CH2-CH-CH3.")
    print("The remaining atoms (NH2) and the phenyl group must be attached to this backbone.")
    print("\nTwo primary isomers are possible:")
    print("  1. 1-phenylpropan-2-amine: C6H5-CH2-CH(NH2)-CH3")
    print("  2. 2-phenylpropan-1-amine: C6H5-CH(CH3)-CH2-NH2")

    # --- Step 4: Isomer Differentiation using 13C NMR ---
    print("\nStep 4: Final Identification")
    c_ch = hsqc_correlations['2.4']
    c_ch2 = hsqc_correlations['2.8']
    print(f"From HSQC, the CH carbon is at {c_ch} ppm and the CH2 carbon is at {c_ch2} ppm.")
    print("In isomer 1, the CH is alpha to the amine (more deshielded) and the CH2 is benzylic.")
    print("In isomer 2, the CH2 is alpha to the amine (more deshielded) and the CH is benzylic.")
    print(f"The data shows C(CH) > C(CH2) ({c_ch} > {c_ch2}), which is consistent with isomer 1.")
    
    # --- Final Answer ---
    print("\n--- Conclusion ---")
    final_name = "1-phenylpropan-2-amine"
    print(f"The compound is identified as {final_name}.")
    print(f"IUPAC Name: {final_name}")

# Execute the analysis
solve_structure()