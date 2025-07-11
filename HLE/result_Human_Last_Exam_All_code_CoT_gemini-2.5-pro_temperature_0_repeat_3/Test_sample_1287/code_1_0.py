def calculate_glycan_mz():
    """
    Calculates the expected m/z for a disialylated, biantennary glycan after
    amidation, permethylation, and sodiation.
    """

    # Monoisotopic atomic masses
    mass_C = 12.000000
    mass_H = 1.007825
    mass_N = 14.003074
    mass_O = 15.994915
    mass_Na = 22.989770

    # Step 1: Calculate the mass of the initial glycan.
    # Composition: 2x GlcNAc, 3x Man, 2x Gal, 2x Neu5Ac. Total 9 residues.
    # Formula for the intact glycan (sum of monosaccharides minus 8x H2O):
    # GlcNAc: C8H15NO6
    # Man/Gal: C6H12O6
    # Neu5Ac: C11H19NO9
    # Total C = (2*8) + (5*6) + (2*11) = 16 + 30 + 22 = 68
    # Total H = (2*15) + (5*12) + (2*19) - (8*2) = 30 + 60 + 38 - 16 = 112
    # Total N = (2*1) + (2*1) = 4
    # Total O = (2*6) + (5*6) + (2*9) - (8*1) = 12 + 30 + 18 - 8 = 52
    # Initial glycan formula: C68 H112 N4 O52
    initial_mass = (68 * mass_C) + (112 * mass_H) + (4 * mass_N) + (52 * mass_O)

    # Step 2: Calculate mass change from amidation.
    # Two sialic acid -COOH groups are converted to -CONH2.
    # Net change per site: -OH + NH2 => -O + NH
    # Mass change = Mass(N) + Mass(H) - Mass(O)
    amidation_change_per_site = mass_N + mass_H - mass_O
    total_amidation_change = 2 * amidation_change_per_site

    # Step 3: Calculate mass change from permethylation.
    # Permethylation replaces -H with -CH3, a net addition of -CH2.
    # We need to count the number of available -OH and -NH sites.
    # Based on the standard A2G2S2 structure:
    # - Reducing GlcNAc: 4 sites (C1,C3,C6-OH, C2-NH)
    # - Core GlcNAc: 3 sites (C3,C6-OH, C2-NH)
    # - Branching Man: 2 sites (C2,C4-OH)
    # - Two outer Man: 3 sites each (6 total)
    # - Two Gal: 3 sites each (6 total)
    # - Two Neu5Ac (amidated): 5 sites each (4 -OH, 1 -NH) (10 total)
    # Total sites = 4 + 3 + 2 + 6 + 6 + 10 = 31 sites.
    # A more robust method: sum(sites_per_monomer - bonds_per_monomer)
    # Monomer sites (all OH/NH): GlcNAc=5, Man/Gal=5, NeuAc=6
    # Bonds: Red.GlcNAc=1, Core.GlcNAc=2, Branch.Man=3, Arm.Man=2, Gal=2, NeuAc=1
    # Sites = (5-1) + (5-2) + (5-3) + 2*(5-2) + 2*(5-2) + 2*(6-1) = 4+3+2+6+6+10 = 31
    num_methylation_sites = 31
    permethylation_change_per_site = mass_C + (2 * mass_H)
    total_permethylation_change = num_methylation_sites * permethylation_change_per_site

    # Step 4: Calculate the final m/z for the [M+Na]+ ion.
    final_neutral_mass = initial_mass + total_amidation_change + total_permethylation_change
    final_mz = final_neutral_mass + mass_Na

    print("The three glycans are isomers and will have the same m/z after derivatization.")
    print("\nCalculation Breakdown:")
    print(f"1. Initial Mass of Glycan (C68H112N4O52): {initial_mass:.4f} Da")
    print(f"2. Mass Change from Amidation (2 sites): {total_amidation_change:.4f} Da")
    print(f"3. Mass Change from Permethylation ({num_methylation_sites} sites): {total_permethylation_change:.4f} Da")
    print(f"4. Mass of Sodium Adduct: {mass_Na:.4f} Da")
    print("\nFinal m/z Calculation:")
    print(f"m/z = {initial_mass:.4f} + ({total_amidation_change:.4f}) + {total_permethylation_change:.4f} + {mass_Na:.4f}")
    print(f"Expected m/z for [M+Na]+ ion: {final_mz:.4f}")

calculate_glycan_mz()
<<<2272.1312>>>