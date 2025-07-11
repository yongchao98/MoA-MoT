def calculate_glycan_mass():
    """
    Calculates the m/z of a permethylated, amidated, singly-sodiated
    A2G2S2 N-glycan.
    """
    # Monoisotopic atomic masses (in Daltons, Da)
    H = 1.007825
    C = 12.000000
    N = 14.003074
    O = 15.994915
    Na = 22.989769

    # Step 1: Base glycan information
    # The native glycan is A2G2S2, which has the formula C84 H138 N6 O62.
    # The three glycans mentioned are isomers and have the same formula and mass.
    print("Step 1: Calculate the mass of the native A2G2S2 glycan (Formula: C84H138N6O62)")
    native_glycan_mass = 84 * C + 138 * H + 6 * N + 62 * O
    print(f"Mass of Native Glycan = (84 * {C}) + (138 * {H}) + (6 * {N}) + (62 * {O}) = {native_glycan_mass:.5f} Da\n")

    # Step 2: Calculate the mass change from amidation
    # Amidation converts -COOH to -CONH2, a net change of -O +N +H.
    # This happens on both of the two sialic acids.
    print("Step 2: Calculate mass change from amidation of two sialic acids")
    amidation_change_per_site = N + H - O
    total_amidation_change = 2 * amidation_change_per_site
    print(f"Mass change per site (-O +N +H) = {N:.5f} + {H:.5f} - {O:.5f} = {amidation_change_per_site:.5f} Da")
    print(f"Total amidation change for 2 sites = 2 * {amidation_change_per_site:.5f} = {total_amidation_change:.5f} Da\n")

    # Step 3: Calculate the mass change from permethylation
    # Permethylation replaces an acidic H with a CH3 group.
    # The net mass change is equivalent to adding a CH2 group.
    # For A2G2S2, there are 37 methylation sites.
    num_methylation_sites = 37
    methylation_change_per_site = C + 2 * H
    total_methylation_change = num_methylation_sites * methylation_change_per_site
    print(f"Step 3: Calculate mass change from permethylation at {num_methylation_sites} sites")
    print(f"Mass change per site (+CH2) = {C:.5f} + (2 * {H:.5f}) = {methylation_change_per_site:.5f} Da")
    print(f"Total methylation change for {num_methylation_sites} sites = {num_methylation_sites} * {methylation_change_per_site:.5f} = {total_methylation_change:.5f} Da\n")

    # Step 4: Calculate the mass of the final derivatized neutral molecule
    final_neutral_mass = native_glycan_mass + total_amidation_change + total_methylation_change
    print("Step 4: Calculate the mass of the neutral derivatized glycan")
    print(f"Final Neutral Mass = Native Mass + Amidation Change + Methylation Change")
    print(f"Final Neutral Mass = {native_glycan_mass:.5f} + ({total_amidation_change:.5f}) + {total_methylation_change:.5f} = {final_neutral_mass:.5f} Da\n")

    # Step 5: Calculate the final m/z of the singly sodiated ion [M+Na]+
    final_sodiated_ion_mz = final_neutral_mass + Na
    print("Step 5: Calculate the m/z of the singly-sodiated ion [M+Na]+")
    print(f"Final m/z = Final Neutral Mass + Mass of Sodium")
    print(f"Final m/z = {final_neutral_mass:.5f} + {Na:.5f} = {final_sodiated_ion_mz:.5f}\n")

    # Final Answer
    print("---")
    print("The calculated m/z is identical for all three isomeric glycans.")
    print(f"The expected mass-to-charge ratio (m/z) for [A2G2S2(amide)2(permethyl)+Na]+ is: {final_sodiated_ion_mz:.4f}")

calculate_glycan_mass()
<<<2762.3840>>>