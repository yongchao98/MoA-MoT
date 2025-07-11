def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized disialylated biantennary N-glycan.
    The glycan undergoes amidation of its sialic acids, followed by permethylation,
    and is detected as a singly sodiated ion.
    """
    # 1. Define monoisotopic masses of atoms and the sodium ion.
    C = 12.000000
    H = 1.007825
    O = 15.994915
    N = 14.003074
    Na = 22.989770

    # 2. Calculate masses of monosaccharide residues and water.
    # Hexose (e.g., Mannose, Galactose): C6H12O6
    mass_hexose = 6 * C + 12 * H + 6 * O
    # N-acetylglucosamine (HexNAc): C8H15NO6
    mass_hexnac = 8 * C + 15 * H + 1 * N + 6 * O
    # N-acetylneuraminic acid (Neu5Ac): C11H19NO9
    mass_neu5ac = 11 * C + 19 * H + 1 * N + 9 * O
    # Water (H2O), lost during glycosidic bond formation
    mass_water = 2 * H + 1 * O

    # 3. Define the glycan composition.
    # A disialylated biantennary glycan (A2G2S2) consists of:
    # 3 Mannose + 2 Galactose = 5 Hexose
    # 4 GlcNAc = 4 HexNAc
    # 2 Neu5Ac = 2 Sialic Acid
    num_hexose = 5
    num_hexnac = 4
    num_neu5ac = 2
    num_residues = num_hexose + num_hexnac + num_neu5ac
    num_linkages = num_residues - 1

    print("--- Step 1: Calculate Initial Mass of Native Glycan ---")
    initial_mass = (num_hexose * mass_hexose) + \
                   (num_hexnac * mass_hexnac) + \
                   (num_neu5ac * mass_neu5ac) - \
                   (num_linkages * mass_water)
    
    print(f"Composition: {num_hexose} Hexose + {num_hexnac} HexNAc + {num_neu5ac} Neu5Ac")
    print(f"Calculation: ({num_hexose} * {mass_hexose:.4f}) + ({num_hexnac} * {mass_hexnac:.4f}) + ({num_neu5ac} * {mass_neu5ac:.4f}) - ({num_linkages} * {mass_water:.4f})")
    print(f"Initial Mass = {initial_mass:.4f} u\n")

    # 4. Calculate mass change from amidation.
    # Reaction: -COOH becomes -CONH2. Net change: -O + NH
    mass_change_amidation = -O + N + H
    num_amidations = num_neu5ac
    mass_after_amidation = initial_mass + (num_amidations * mass_change_amidation)
    
    print("--- Step 2: Calculate Mass After Amidation ---")
    print(f"Amidation converts {num_amidations} COOH groups to CONH2 groups.")
    print(f"Mass change per amidation: {mass_change_amidation:.4f} u")
    print(f"Calculation: {initial_mass:.4f} + ({num_amidations} * {mass_change_amidation:.4f})")
    print(f"Mass after Amidation = {mass_after_amidation:.4f} u\n")

    # 5. Calculate the number of permethylation sites.
    # Permethylation adds a methyl group to each -OH and -NH group.
    # Sites per monomer: Hexose=5, HexNAc=4 (3 OH, 1 NH), Amidated Neu5Ac=8 (5 OH, 1 NH, 1 CONH2)
    sites_hex = 5
    sites_hexnac = 4
    sites_neu5ac_amide = 8
    total_potential_sites = (num_hexose * sites_hex) + (num_hexnac * sites_hexnac) + (num_neu5ac * sites_neu5ac_amide)
    sites_lost_to_linkages = 2 * num_linkages
    num_methylation_sites = total_potential_sites - sites_lost_to_linkages
    
    print("--- Step 3: Calculate Number of Permethylation Sites ---")
    print(f"Total potential sites: ({num_hexose} * {sites_hex}) + ({num_hexnac} * {sites_hexnac}) + ({num_neu5ac} * {sites_neu5ac_amide}) = {total_potential_sites}")
    print(f"Sites lost to {num_linkages} linkages: 2 * {num_linkages} = {sites_lost_to_linkages}")
    print(f"Calculation: {total_potential_sites} - {sites_lost_to_linkages}")
    print(f"Total Permethylation Sites = {num_methylation_sites}\n")

    # 6. Calculate mass change from permethylation.
    # Reaction: -H becomes -CH3. Net change: +CH2
    mass_change_methylation = C + 2 * H
    mass_after_permethylation = mass_after_amidation + (num_methylation_sites * mass_change_methylation)

    print("--- Step 4: Calculate Mass After Permethylation ---")
    print(f"Permethylation adds {num_methylation_sites} methyl groups.")
    print(f"Mass increase per site: {mass_change_methylation:.4f} u")
    print(f"Calculation: {mass_after_amidation:.4f} + ({num_methylation_sites} * {mass_change_methylation:.4f})")
    print(f"Final Neutral Mass = {mass_after_permethylation:.4f} u\n")

    # 7. Calculate final m/z of the sodiated ion [M+Na]+.
    final_mz = mass_after_permethylation + Na
    
    print("--- Step 5: Calculate Final m/z for Singly Sodiated Ion [M+Na]+ ---")
    print(f"Calculation: {mass_after_permethylation:.4f} (M) + {Na:.4f} (Na+)")
    print(f"Final m/z [M+Na]+ = {final_mz:.4f}\n")
    
    print("--- Final Answer ---")
    print("The three glycans are isomers, and the chemical modifications are not linkage-specific.")
    print("Therefore, all three will have the same expected mass.")
    print(f"Mass for A2G(4)2S(3)2: {final_mz:.4f}")
    print(f"Mass for A2G(4)S(3)S(6): {final_mz:.4f}")
    print(f"Mass for A2G(4)2S(6)2: {final_mz:.4f}")
    
    return final_mz

final_mass_to_charge = calculate_glycan_mass()
# The final answer is wrapped in <<<>>> as requested.
# The value is formatted to 4 decimal places as is common in mass spectrometry.
print(f"\n<<<{final_mass_to_charge:.4f}>>>")
