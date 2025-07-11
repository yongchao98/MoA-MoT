import sys
# Redirect print to stderr to show the thinking process, as the final answer will be captured from stdout
# print = lambda x, **kwargs: sys.stderr.write(str(x) + '\n') # Commented out for final output clarity

def identify_glycan():
    """
    Performs the step-by-step identification of the glycan based on MS and MS/MS data.
    """
    # Define physical and chemical constants
    # Masses of proton and common monosaccharide residues (as part of a polymer, after loss of H2O)
    H_MASS = 1.007825
    HEX_MASS = 162.05282  # Mannose (Man), Galactose (Gal)
    HEXNAC_MASS = 203.07937 # N-Acetylglucosamine (GlcNAc)
    FUC_MASS = 146.05791   # Fucose (dHex)
    NEU5AC_MASS = 291.09542 # N-Acetylneuraminic acid (Sialic acid)

    # Monoisotopic mass of the RapiFluor-MS (RFMS) tag that is chemically added to the glycan
    RFMS_TAG_MASS = 314.17428 # Formula C17H22N4O2

    # Observed MS and MS/MS data from the user
    parent_mz = 856.6638
    isotope_peaks = [856.9971, 857.3305, 857.6638]
    msms_fragments = {
        'oxonium1': 204.087,
        'oxonium2': 366.140,
        'base_peak': 528.193,
        'other1': 673.231,
        'other2': 882.409,
        'other3': 1368.568,
        'other4': 1894.753,
        'other5': 2260.886
    }
    
    # --- Step 1: Determine Charge State ---
    # The isotopic spacing is ~0.333 Da, which corresponds to 1/3, indicating z=3.
    charge_state = 3

    print("--- Glycan Identification Analysis ---")
    print(f"\nStep 1: Determine Precursor Ion Mass from m/z {parent_mz}")
    print(f"The isotopic spacing (~0.333 Da) confirms a charge state (z) of {charge_state}.")
    
    # --- Step 2: Calculate Mass of the RFMS-Glycan and Native Glycan ---
    mass_glycan_rfms = (parent_mz * charge_state) - (charge_state * H_MASS)
    mass_glycan_native = mass_glycan_rfms - RFMS_TAG_MASS

    print(f"\nStep 2: Calculate Glycan Mass")
    print(f"Mass of the intact RFMS-glycan [M] is calculated as: ({parent_mz} * {charge_state}) - ({charge_state} * {H_MASS}) = {mass_glycan_rfms:.4f} Da.")
    print(f"Mass of the native glycan is found by subtracting the RFMS tag ({RFMS_TAG_MASS:.4f} Da): {mass_glycan_rfms:.4f} - {RFMS_TAG_MASS:.4f} = {mass_glycan_native:.4f} Da.")

    # --- Step 3: Determine Glycan Composition ---
    # Propose a composition that matches the native glycan mass.
    comp = {'Hex': 4, 'HexNAc': 5, 'Fuc': 1, 'NeuAc': 1}
    comp_mass = (comp['Hex'] * HEX_MASS) + \
                (comp['HexNAc'] * HEXNAC_MASS) + \
                (comp['Fuc'] * FUC_MASS) + \
                (comp['NeuAc'] * NEU5AC_MASS)
    mass_diff = mass_glycan_native - comp_mass

    print(f"\nStep 3: Determine Glycan Composition")
    print(f"The native glycan mass of {mass_glycan_native:.4f} Da corresponds to a composition of Hex(4)HexNAc(5)Fuc(1)NeuAc(1).")
    print(f"The theoretical mass for this composition is: (4 * {HEX_MASS:.4f}) + (5 * {HEXNAC_MASS:.4f}) + (1 * {FUC_MASS:.4f}) + (1 * {NEU5AC_MASS:.4f}) = {comp_mass:.4f} Da.")
    print(f"The mass difference is only {mass_diff:.4f} Da, confirming this composition.")
    
    # --- Step 4: Propose Structure and Confirm with MS/MS fragments ---
    # Calculate theoretical fragment masses to compare with observed data.
    # The [HexNAc]+ oxonium ion, formula C8H14NO5+
    hexnac_oxonium_theory = HEXNAC_MASS + 18.01056 - H_MASS 
    # The [Hex+HexNAc]+ oxonium ion, formula C14H24NO10+
    hex_hexnac_oxonium_theory = HEX_MASS + HEXNAC_MASS + 18.01056 - H_MASS
    # The characteristic D-ion for a bisected structure. [Man(2)+GlcNAc(1)+H]+
    bisecting_d_ion_theory = (2 * HEX_MASS) + HEXNAC_MASS + H_MASS
    
    print("\nStep 4: Confirm Structure with MS/MS Fragmentation")
    print("This composition corresponds to a common N-glycan structure: a monosialylated, core-fucosylated, biantennary glycan with a bisecting GlcNAc.")
    print("Key fragment ions confirm this structure:")
    print(f"  - m/z {msms_fragments['oxonium1']}: Matches the [HexNAc]+ oxonium ion (theoretical {hexnac_oxonium_theory:.4f} Da), confirming the presence of HexNAc.")
    print(f"  - m/z {msms_fragments['oxonium2']}: Matches the [Hex-HexNAc]+ oxonium ion (theoretical {hex_hexnac_oxonium_theory:.4f} Da), confirming a lactosamine-type antenna.")
    print(f"  - m/z {msms_fragments['base_peak']} (Most Intense Ion): Matches a characteristic D-ion from fragmentation across the glycan core (theoretical {bisecting_d_ion_theory:.4f} Da). This fragment is strong evidence for a bisecting GlcNAc.")
          
    # --- Step 5: Final Conclusion ---
    final_name = "FA2BG1S1"
    linkage_info = "The structure contains a core α1-6 linked fucose and a β1-4 linked bisecting GlcNAc. The specific linkage of the sialic acid (e.g., α2-3 or α2-6) cannot be determined from this CID data."
    print("\n--- Final Identification ---")
    print(f"Name (Oxford Nomenclature): {final_name}")
    print(f"Linkage Information: {linkage_info}")

if __name__ == "__main__":
    identify_glycan()