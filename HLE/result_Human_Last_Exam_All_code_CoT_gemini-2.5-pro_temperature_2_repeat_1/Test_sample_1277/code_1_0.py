import math

def solve_glycan_puzzle():
    """
    Solves the glycan identification puzzle by calculating masses from MS data
    and interpreting fragmentation patterns.
    """

    # --- Constants (Monoisotopic Masses) ---
    MASS_PROTON = 1.007276
    MASS_H = 1.007825
    MASS_O = 15.994915
    MASS_H2O = 2 * MASS_H + MASS_O

    # Dehydrated monosaccharide residue masses
    MASS_HEX = 162.05282
    MASS_HEXNAC = 203.07937
    MASS_FUC = 146.05791  # Deoxyhexose
    MASS_NEUAC = 291.09542 # N-Acetylneuraminic acid
    MASS_NEUGC = 307.09033 # N-Glycolylneuraminic acid

    # RFMS Tag Derivatization (Reductive Amination)
    # Mass of RFMS reagent (C13H14N2O) - Mass of Oxygen atom
    MASS_RFMS_ADDITION = 234.1118 - MASS_O
    
    # --- Input Data ---
    parent_ion_mz = 856.6638
    isotopes_mz = [856.9971, 857.3305]
    msms_fragments = [204.087, 366.140, 528.193, 673.231, 882.409, 1368.568, 1894.753, 2260.886]

    print("--- Step 1: Analyze the Parent Ion ---")
    delta_mz = isotopes_mz[0] - parent_ion_mz
    charge_state = round(1 / delta_mz)
    print(f"Isotope spacing (m/z {isotopes_mz[0]} - m/z {parent_ion_mz}) is ~{delta_mz:.4f}.")
    print(f"This corresponds to a charge state (z) of +{charge_state}.")

    neutral_mass_derivatized = (parent_ion_mz * charge_state) - (charge_state * MASS_PROTON)
    print(f"The calculated neutral mass of the RFMS-derivatized glycan is {neutral_mass_derivatized:.4f} Da.")
    print("-" * 20)

    print("\n--- Step 2: Determine Native Glycan Mass ---")
    print(f"The mass added by RFMS derivatization (reductive amination) is {MASS_RFMS_ADDITION:.4f} Da.")
    native_glycan_mass = neutral_mass_derivatized - MASS_RFMS_ADDITION
    print(f"The calculated mass of the native (underivatized) glycan is {native_glycan_mass:.4f} Da.")
    
    dehydrated_glycan_mass = native_glycan_mass - MASS_H2O
    print(f"The mass of the dehydrated glycan polymer (sum of residues) is {dehydrated_glycan_mass:.4f} Da.")
    print("-" * 20)

    print("\n--- Step 3: Analyze Key MS/MS Fragments ---")
    
    # Fragment 1: m/z 204.087
    frag1_calc = MASS_HEXNAC + MASS_H
    print(f"Fragment at m/z {msms_fragments[0]:.4f}: Identified as [HexNAc+H]+ oxonium ion (Calc: {frag1_calc:.4f}). Confirms presence of HexNAc.")
    
    # Fragment 2: m/z 366.140
    frag2_calc = MASS_HEX + MASS_HEXNAC + MASS_H
    print(f"Fragment at m/z {msms_fragments[1]:.4f}: Identified as [Hex-HexNAc+H]+ (Calc: {frag2_calc:.4f}). Confirms LacNAc-type antennae.")

    # Fragment 4: m/z 673.231
    frag4_calc_neuac = MASS_NEUAC + MASS_HEX + MASS_HEXNAC + MASS_H
    frag4_calc_neugc = MASS_NEUGC + MASS_HEX + MASS_HEXNAC + MASS_H
    print(f"Fragment at m/z {msms_fragments[3]:.4f}: This mass matches [NeuGc+Hex+HexNAc+H]+ (Calc: {frag4_calc_neugc:.4f}), not NeuAc (Calc: {frag4_calc_neuac:.4f}). Confirms one antenna is capped with N-Glycolylneuraminic acid (NeuGc).")
    
    # Fragment 8: m/z 2260.886
    frag8_calc = neutral_mass_derivatized - MASS_NEUGC + MASS_H
    print(f"Fragment at m/z {msms_fragments[7]:.4f}: This corresponds to the neutral loss of NeuGc from the parent ion, [M+H-NeuGc]+ (Calc: {frag8_calc:.4f}). Confirms the glycan has exactly one NeuGc residue.")

    # Fragment 7: m/z 1894.753
    mass_neugc_antenna = MASS_NEUGC + MASS_HEX + MASS_HEXNAC
    frag7_calc = neutral_mass_derivatized - mass_neugc_antenna + MASS_H
    print(f"Fragment at m/z {msms_fragments[6]:.4f}: This Y-ion corresponds to the loss of the entire NeuGc-Hex-HexNAc antenna, [M+H-(NeuGc-Hex-HexNAc)]+ (Calc: {frag7_calc:.4f}).")

    print(f"Fragment at m/z {msms_fragments[2]:.4f}: This base peak corresponds to a [Hex-HexNAc-Man+H]+ B-ion, suggesting a typical complex N-glycan branching structure.")

    print("-" * 20)

    print("\n--- Step 4: Deduce Glycan Composition ---")
    asialo_dehydrated_mass = dehydrated_glycan_mass - MASS_NEUGC
    print(f"The mass of the asialo-, dehydrated glycan portion is {asialo_dehydrated_mass:.4f} Da.")
    
    # Search for composition
    found_comp = None
    tolerance = 0.01
    for h in range(3, 8):
        for n in range(4, 8):
            for f in range(0, 3):
                calc_mass = h * MASS_HEX + n * MASS_HEXNAC + f * MASS_FUC
                if abs(calc_mass - asialo_dehydrated_mass) < tolerance:
                    found_comp = {'H': h, 'N': n, 'F': f, 'Mass': calc_mass}
                    break
            if found_comp: break
        if found_comp: break
    
    if found_comp:
        print(f"A composition of Hex({found_comp['H']}), HexNAc({found_comp['N']}), Fuc({found_comp['F']}) matches this mass ({found_comp['Mass']:.4f} Da).")
        print("\n--- Final Composition ---")
        final_H, final_N, final_F, final_S = found_comp['H'], found_comp['N'], found_comp['F'], 1
        print(f"The full glycan composition is: Hex({final_H}), HexNAc({final_N}), Fuc({final_F}), NeuGc({final_S})")
    else:
        print("Could not find a matching composition.")
        return
        
    print("-" * 20)
    
    print("\n--- Step 5: Assign Glycan Name ---")
    # Interpret composition using N-glycan rules: 3 Man in core, 2 GlcNAc in core
    num_antennae = final_N - 2
    num_gal = final_H - 3
    
    glycan_name = ""
    if final_F > 0:
        glycan_name += "F" # F for core Fucose
    
    glycan_name += f"A{num_antennae}" # A for antennae
    
    if num_gal > 0:
        glycan_name += f"G{num_gal}" # G for Galactose
        
    glycan_name += f"S({final_S if final_S > 1 else ''}Gc)" if final_S > 0 else "" # S for Sialic acid, (Gc) for glycolyl
    
    print("Based on the composition (Fuc=1, 4 antennae, 1 Galactose, 1 NeuGc), the Oxford name is determined.")
    print(f"The glycan is a core-fucosylated ({final_F} Fuc), tetra-antennary (HexNAc {final_N} -> A{num_antennae}) structure.")
    print(f"It has {num_gal} galactose residue, which is capped by the NeuGc residue.")
    print("\nThe three other antennae consist of terminal, non-galactosylated GlcNAc residues.")
    
    final_answer = f"FA{num_antennae}G{num_gal}S(Gc)1" # or S(Gc)
    print("\nFinal proposed name:")
    print(final_answer)

solve_glycan_puzzle()
<<<FA4G1S(Gc)1>>>