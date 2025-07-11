import sys
# Helper function to prevent python from scientific notation for cleaner output
def format_float(num):
    return f'{num:.4f}'

def calculate_glycan_mass():
    """
    Calculates the expected m/z for specific derivatized N-glycans.

    The process involves three steps:
    1. Amidation: Sialic acid -COOH groups are converted to -CONH2.
    2. Permethylation: All -OH and -NH hydrogens are replaced with methyl groups.
    3. Sodiated Ion Formation: A sodium ion adducts for MS detection.

    The three glycans mentioned (A2G(4)2S(3)2, A2G(4)S(3)S(6), A2G(4)2S(6)2)
    are isomers and will have the identical mass after these reactions.
    """
    # Monoisotopic atomic masses
    ATOM_MASS = {
        'C': 12.0000000,
        'H': 1.0078250,
        'O': 15.9949146,
        'N': 14.0030740,
        'Na': 22.989770,
    }

    # Calculate masses of fundamental units
    mass_h2o = 2 * ATOM_MASS['H'] + ATOM_MASS['O']
    mass_hexose = 6 * ATOM_MASS['C'] + 12 * ATOM_MASS['H'] + 6 * ATOM_MASS['O']
    mass_hexnac = 8 * ATOM_MASS['C'] + 15 * ATOM_MASS['H'] + ATOM_MASS['N'] + 6 * ATOM_MASS['O']
    mass_neuac = 11 * ATOM_MASS['C'] + 19 * ATOM_MASS['H'] + ATOM_MASS['N'] + 9 * ATOM_MASS['O']

    # --- Step 1: Calculate the mass of the native A2G2S2 glycan ---
    # Composition: 5 Hexose, 4 HexNAc, 2 NeuAc. Total 11 residues, 10 glycosidic bonds.
    # The mass is the sum of monomers minus the mass of water for each bond formed.
    num_hex = 5
    num_hexnac = 4
    num_neuac = 2
    num_bonds = num_hex + num_hexnac + num_neuac - 1
    
    mass_native_glycan = (num_hex * mass_hexose +
                          num_hexnac * mass_hexnac +
                          num_neuac * mass_neuac) - num_bonds * mass_h2o

    # --- Step 2: Calculate mass change from Amidation ---
    # Reaction on carboxyl group: -COOH -> -CONH2.
    # This is equivalent to replacing one -OH group with an -NH2 group.
    mass_oh = ATOM_MASS['O'] + ATOM_MASS['H']
    mass_nh2 = ATOM_MASS['N'] + 2 * ATOM_MASS['H']
    delta_amidation_per_site = mass_nh2 - mass_oh
    
    num_sialic_acids = 2
    total_delta_amidation = num_sialic_acids * delta_amidation_per_site
    
    mass_amidated_glycan = mass_native_glycan + total_delta_amidation

    # --- Step 3: Calculate mass change from Permethylation ---
    # Determine the number of methylation sites on the amidated glycan.
    # We start with the known number of sites for the A2G2 core (37)
    # and add the sites from the two attached, amidated sialic acid residues.
    sites_A2G2_core = 37
    # Each amidated sialic acid residue has 4 -OH, 1 N-acetyl -NH, and 1 amide -NH2 -> 4+1+2=7 sites.
    # Adding it to the core uses up one -OH site on a galactose.
    # Net sites added per residue = 7 (from Sia-amide) - 1 (from Gal) = 6.
    net_sites_added = (7 - 1) * num_sialic_acids
    total_methylation_sites = sites_A2G2_core + net_sites_added
    
    # Mass change per methylation: replacing -H with -CH3, a net addition of -CH2.
    delta_methylation_per_site = ATOM_MASS['C'] + 2 * ATOM_MASS['H']
    total_delta_methylation = total_methylation_sites * delta_methylation_per_site
    
    mass_derivatized_glycan = mass_amidated_glycan + total_delta_methylation
    
    # --- Step 4: Calculate final m/z of the sodiated ion [M+Na]+ ---
    mass_sodium = ATOM_MASS['Na']
    final_mz = mass_derivatized_glycan + mass_sodium

    # --- Output the results ---
    print("Since A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers, they all have the same mass.")
    print("The expected m/z value is identical for all three compounds.\n")
    print("Calculation Steps:")
    print(f"1. Initial Mass of A2G2S2 Glycan [M]: {format_float(mass_native_glycan)} Da")
    print(f"2. Mass after Amidation [M_amidated] = {format_float(mass_native_glycan)} + ({num_sialic_acids} * {format_float(delta_amidation_per_site)}) = {format_float(mass_amidated_glycan)} Da")
    print(f"3. Mass after Permethylation [M'] = {format_float(mass_amidated_glycan)} + ({total_methylation_sites} sites * {format_float(delta_methylation_per_site)}) = {format_float(mass_derivatized_glycan)} Da")
    print("-" * 30)
    print("Final Mass Calculation:")
    print(f"m/z of [M'+Na]+ = {format_float(mass_derivatized_glycan)} (M') + {format_float(mass_sodium)} (Na+)")
    print("-" * 30)
    
    final_answer_text = (
        f"The expected mass for A2G(4)2S(3)2 is {format_float(final_mz)} m/z.\n"
        f"The expected mass for A2G(4)S(3)S(6) is {format_float(final_mz)} m/z.\n"
        f"The expected mass for A2G(4)2S(6)2 is {format_float(final_mz)} m/z."
    )
    print(final_answer_text)

    # Outputting the final answer in the specified format
    # Redirecting to stderr to not interfere with stdout for a clean final number.
    print(f'<<<{format_float(final_mz)}>>>', file=sys.stderr)

calculate_glycan_mass()
