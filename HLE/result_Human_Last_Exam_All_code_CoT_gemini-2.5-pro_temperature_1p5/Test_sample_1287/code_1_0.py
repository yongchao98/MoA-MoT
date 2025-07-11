def calculate_glycan_masses():
    """
    Calculates the theoretical m/z values for three derivatized sialylated glycans.

    The calculation is based on the hypothesis that under the specified reaction
    conditions, α2,6-linked sialic acids are amidated, while α2,3-linked
    sialic acids undergo lactonization, prior to permethylation.
    """

    # Monoisotopic atomic masses
    H = 1.007825
    C = 12.000000
    N = 14.003074
    O = 15.994915
    Na = 22.98977

    # --- Step 1: Base Glycan Composition and Mass ---
    # Composition: 5 Hexose (Man/Gal), 4 GlcNAc, 2 Neu5Ac
    # Formula: C84 H138 N6 O62
    base_glycan_formula = {'C': 84, 'H': 138, 'N': 6, 'O': 62}
    base_glycan_mass = (base_glycan_formula['C'] * C +
                        base_glycan_formula['H'] * H +
                        base_glycan_formula['N'] * N +
                        base_glycan_formula['O'] * O)

    # --- Step 2: Mass Changes from Reactions ---
    # Mass gain per methylation (replaces H with CH3 -> adds CH2)
    permethylation_gain = C + 2 * H # 14.01565 Da
    
    # Mass change for amidation (-COOH -> -CONH2; effectively -OH + NH2)
    amidation_change = (N + 2 * H) - O # -0.984016 Da
    
    # Mass change for lactonization (loss of H2O)
    lactonization_change = - (2 * H + O) # -18.010565 Da

    # --- Step 3: Count Methylation Sites ---
    # For a free reducing A2G2S2 N-glycan, there are 31 -OH groups and 6 N-acetyl -NH groups.
    # Total "backbone" sites = 37.
    backbone_sites = 37
    # Amidation of -COOH to -CONH2 adds 2 methylation sites (on the -NH2).
    amide_sites_added = 2
    # Lactonization consumes one of the sialic acid's own -OH groups.
    lactone_sites_lost = 1

    print("Analysis of Derivatized Sialylated Glycans\n")
    print(f"Starting Glycan (A2G2S2) Mass: {base_glycan_mass:.4f} Da")
    print("-" * 50)

    # --- Calculation for A2G(4)2S(6)2 (Amide x 2) ---
    glycan_1_name = "A2G(4)2S(6)2"
    num_amides_1 = 2
    num_lactones_1 = 0
    reaction_change_1 = num_amides_1 * amidation_change
    methyl_sites_1 = backbone_sites + num_amides_1 * amide_sites_added
    methylation_mass_1 = methyl_sites_1 * permethylation_gain
    final_mass_1 = base_glycan_mass + reaction_change_1 + methylation_mass_1
    final_mz_1 = final_mass_1 + Na

    print(f"For {glycan_1_name} (assumed Amide x2):")
    print(f"m/z = (Base Mass + 2 * Amidation Change + {methyl_sites_1} * Permethylation Gain) + Na Mass")
    print(f"m/z = ({base_glycan_mass:.4f} + 2 * {amidation_change:.4f} + {methyl_sites_1} * {permethylation_gain:.4f}) + {Na:.4f}")
    print(f"m/z = ({base_glycan_mass:.4f} {reaction_change_1:+.4f} {methylation_mass_1:+.4f}) + {Na:.4f}")
    print(f"Expected m/z for [M+Na]+: {final_mz_1:.4f}\n")

    # --- Calculation for A2G(4)S(3)S(6) (Amide x 1, Lactone x 1) ---
    glycan_2_name = "A2G(4)S(3)S(6)"
    num_amides_2 = 1
    num_lactones_2 = 1
    reaction_change_2 = (num_amides_2 * amidation_change) + (num_lactones_2 * lactonization_change)
    methyl_sites_2 = backbone_sites + (num_amides_2 * amide_sites_added) - (num_lactones_2 * lactone_sites_lost)
    methylation_mass_2 = methyl_sites_2 * permethylation_gain
    final_mass_2 = base_glycan_mass + reaction_change_2 + methylation_mass_2
    final_mz_2 = final_mass_2 + Na
    
    print(f"For {glycan_2_name} (assumed Amide x1, Lactone x1):")
    print(f"m/z = (Base Mass + Amidation Change + Lactonization Change + {methyl_sites_2} * Permethylation Gain) + Na Mass")
    print(f"m/z = ({base_glycan_mass:.4f} + {amidation_change:.4f} + {lactonization_change:.4f} + {methyl_sites_2} * {permethylation_gain:.4f}) + {Na:.4f}")
    print(f"m/z = ({base_glycan_mass:.4f} {reaction_change_2:+.4f} {methylation_mass_2:+.4f}) + {Na:.4f}")
    print(f"Expected m/z for [M+Na]+: {final_mz_2:.4f}\n")

    # --- Calculation for A2G(4)2S(3)2 (Lactone x 2) ---
    glycan_3_name = "A2G(4)2S(3)2"
    num_amides_3 = 0
    num_lactones_3 = 2
    reaction_change_3 = num_lactones_3 * lactonization_change
    methyl_sites_3 = backbone_sites - (num_lactones_3 * lactone_sites_lost)
    methylation_mass_3 = methyl_sites_3 * permethylation_gain
    final_mass_3 = base_glycan_mass + reaction_change_3 + methylation_mass_3
    final_mz_3 = final_mass_3 + Na

    print(f"For {glycan_3_name} (assumed Lactone x2):")
    print(f"m/z = (Base Mass + 2 * Lactonization Change + {methyl_sites_3} * Permethylation Gain) + Na Mass")
    print(f"m/z = ({base_glycan_mass:.4f} + 2 * {lactonization_change:.4f} + {methyl_sites_3} * {permethylation_gain:.4f}) + {Na:.4f}")
    print(f"m/z = ({base_glycan_mass:.4f} {reaction_change_3:+.4f} {methylation_mass_3:+.4f}) + {Na:.4f}")
    print(f"Expected m/z for [M+Na]+: {final_mz_3:.4f}\n")

# Execute the calculation
calculate_glycan_masses()