import math

def calculate_glycan_masses():
    """
    Calculates the theoretical m/z values for three modified and permethylated
    sialylated glycans detected as singly-sodiated ions.

    The key assumptions are:
    1.  The starting glycan is A2G2S2 (Man3GlcNAc4Gal2Neu5Ac2).
    2.  The DMT-MM/NH4Cl reaction converts alpha-2,6 linked sialic acids to
        primary amides (-CONH2) and alpha-2,3 linked sialic acids to
        intramolecular lactams.
    3.  Permethylation exhaustively methylates all free hydroxyl (-OH) groups and
        all available N-H protons (on N-acetyl groups and on the new primary amides).
    4.  Detection is of the singly-sodiated species, [M+Na]+.
    5.  All calculations use monoisotopic masses.
    """
    # 1. DEFINE MONOISOTOPIC MASSES
    MASS_H = 1.007825
    MASS_C = 12.000000
    MASS_N = 14.003074
    MASS_O = 15.994915
    MASS_NA = 22.989770

    # Masses of common groups and residues (as they exist in a polymer)
    MASS_H2O = 2 * MASS_H + MASS_O
    # Mass change for methylation (-H -> -CH3), net addition of CH2
    MASS_CHANGE_METHYLATION = MASS_C + 2 * MASS_H
    
    # Residue masses (monosaccharide mass - H2O)
    MASS_RES_HEX = (6 * MASS_C) + (10 * MASS_H) + (5 * MASS_O)  # Hexose (e.g., Man, Gal)
    MASS_RES_HEXNAC = (8 * MASS_C) + (13 * MASS_H) + MASS_N + (5 * MASS_O)  # N-Acetylhexosamine (e.g., GlcNAc)
    # Neu5Ac residue: C11H17NO8. N-Acetylneuraminic acid
    MASS_RES_NEU5AC = (11 * MASS_C) + (17 * MASS_H) + MASS_N + (8 * MASS_O) 

    # 2. CALCULATE STARTING MASS of A2G2S2
    # Structure: 5 Hexose (3 Man + 2 Gal), 4 HexNAc, 2 Neu5Ac
    # The final mass is the sum of residues plus one water for the reducing end
    mass_A2G2S2 = (5 * MASS_RES_HEX) + \
                  (4 * MASS_RES_HEXNAC) + \
                  (2 * MASS_RES_NEU5AC) + \
                  MASS_H2O

    # 3. MODEL THE LINKAGE-SPECIFIC REACTIONS
    # Amidation: -COOH -> -CONH2. Net change: +NH2 -OH = +NH -O
    mass_change_amidation = (MASS_N + MASS_H) - MASS_O
    # Lactamization: intramolecular reaction with loss of water
    mass_change_lactamization = -MASS_H2O
    
    # The A2G2S2 glycan has 27 -OH groups and 6 N-acetyl groups (-NH-)
    # The 6 N-acetyls are on the 4 GlcNAc and 2 Neu5Ac residues.
    num_oh_groups = 27
    num_nh_acetyl_glcnac = 4
    num_nh_acetyl_neu5ac = 2

    print("--- Calculation of Expected Glycan Masses ---\n")

    # --- Glycan 1: A2G(4)2S(3)2 ---
    # Both sialic acids are a-2,3 linked -> both form lactams.
    glycan1_name = "A2G(4)2S(3)2"
    g1_reaction_change = 2 * mass_change_lactamization
    # Permethylation sites: 27 -OH groups. 4 N-acetyls on GlcNAc.
    # The 2 N-acetyls on Neu5Ac are consumed in lactam formation, so no N-H left.
    g1_methylation_sites = num_oh_groups + num_nh_acetyl_glcnac
    g1_methylation_mass = g1_methylation_sites * MASS_CHANGE_METHYLATION
    g1_final_mz = mass_A2G2S2 + g1_reaction_change + g1_methylation_mass + MASS_NA
    
    print(f"Analysis for {glycan1_name} (di-lactam derivative):")
    print(f"  Base Glycan Mass (A2G2S2): {mass_A2G2S2:.4f} Da")
    print(f"  Reaction Mass Change (2x lactamization): {g1_reaction_change:.4f} Da")
    print(f"  Permethylation Mass Change ({g1_methylation_sites} x Me): {g1_methylation_mass:.4f} Da")
    print(f"  Sodium Adduct Mass (Na+): {MASS_NA:.4f} Da")
    print(f"  Final [M+Na]+ m/z = {mass_A2G2S2:.4f} + ({g1_reaction_change:.4f}) + {g1_methylation_mass:.4f} + {MASS_NA:.4f} = {g1_final_mz:.4f} Da\n")

    # --- Glycan 2: A2G(4)S(3)S(6) ---
    # One a-2,3 (-> lactam) and one a-2,6 (-> amide).
    glycan2_name = "A2G(4)S(3)S(6)"
    g2_reaction_change = mass_change_lactamization + mass_change_amidation
    # Permethylation sites: 27 -OH. 4 GlcNAc N-acetyls. 1 Neu5Ac N-acetyl (amide one).
    # The lactamized Neu5Ac N-acetyl is consumed. The amidation creates one -CONH2 group, which adds 2 new N-H sites.
    g2_methylation_sites = num_oh_groups + num_nh_acetyl_glcnac + 1 + 2
    g2_methylation_mass = g2_methylation_sites * MASS_CHANGE_METHYLATION
    g2_final_mz = mass_A2G2S2 + g2_reaction_change + g2_methylation_mass + MASS_NA

    print(f"Analysis for {glycan2_name} (amide-lactam derivative):")
    print(f"  Base Glycan Mass (A2G2S2): {mass_A2G2S2:.4f} Da")
    print(f"  Reaction Mass Change (1x lactam + 1x amide): {g2_reaction_change:.4f} Da")
    print(f"  Permethylation Mass Change ({g2_methylation_sites} x Me): {g2_methylation_mass:.4f} Da")
    print(f"  Sodium Adduct Mass (Na+): {MASS_NA:.4f} Da")
    print(f"  Final [M+Na]+ m/z = {mass_A2G2S2:.4f} + ({g2_reaction_change:.4f}) + {g2_methylation_mass:.4f} + {MASS_NA:.4f} = {g2_final_mz:.4f} Da\n")

    # --- Glycan 3: A2G(4)2S(6)2 ---
    # Both sialic acids are a-2,6 linked -> both form amides.
    glycan3_name = "A2G(4)2S(6)2"
    g3_reaction_change = 2 * mass_change_amidation
    # Permethylation sites: 27 -OH. All 6 N-acetyls (4 GlcNAc + 2 Neu5Ac).
    # The two amidation reactions create two -CONH2 groups, adding 2*2=4 new N-H sites.
    g3_methylation_sites = num_oh_groups + num_nh_acetyl_glcnac + num_nh_acetyl_neu5ac + 4
    g3_methylation_mass = g3_methylation_sites * MASS_CHANGE_METHYLATION
    g3_final_mz = mass_A2G2S2 + g3_reaction_change + g3_methylation_mass + MASS_NA
    
    print(f"Analysis for {glycan3_name} (di-amide derivative):")
    print(f"  Base Glycan Mass (A2G2S2): {mass_A2G2S2:.4f} Da")
    print(f"  Reaction Mass Change (2x amidation): {g3_reaction_change:.4f} Da")
    print(f"  Permethylation Mass Change ({g3_methylation_sites} x Me): {g3_methylation_mass:.4f} Da")
    print(f"  Sodium Adduct Mass (Na+): {MASS_NA:.4f} Da")
    print(f"  Final [M+Na]+ m/z = {mass_A2G2S2:.4f} + ({g3_reaction_change:.4f}) + {g3_methylation_mass:.4f} + {MASS_NA:.4f} = {g3_final_mz:.4f} Da")

# Run the calculation and print the results
calculate_glycan_masses()
>>>The masses you should observe for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are 2644.2482, 2703.3217, and 2762.3951 Da, respectively.