import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized glycan based on its composition and modifications.
    """
    # Monoisotopic atomic masses
    ATOMIC_MASS = {
        'H': 1.007825,
        'C': 12.000000,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Monoisotopic masses of monosaccharide residues (as C_x H_y O_z, not including losses for bonds)
    # Galactose is an isomer of Mannose, so they have the same mass.
    MASS_HEX = 6 * ATOMIC_MASS['C'] + 12 * ATOMIC_MASS['H'] + 6 * ATOMIC_MASS['O']
    MASS_HEXNAC = 8 * ATOMIC_MASS['C'] + 15 * ATOMIC_MASS['H'] + 1 * ATOMIC_MASS['N'] + 6 * ATOMIC_MASS['O']
    MASS_NEU5AC = 11 * ATOMIC_MASS['C'] + 19 * ATOMIC_MASS['H'] + 1 * ATOMIC_MASS['N'] + 9 * ATOMIC_MASS['O']
    MASS_H2O = 2 * ATOMIC_MASS['H'] + ATOMIC_MASS['O']

    # Glycan composition
    composition = {
        'Hex': 5,      # 3 Man + 2 Gal
        'HexNAc': 4,   # 4 GlcNAc
        'Neu5Ac': 2,   # 2 Sialic Acid
    }
    
    # --- Step 1: Calculate the mass of the native glycan ---
    num_monosaccharides = sum(composition.values())
    num_glycosidic_bonds = num_monosaccharides - 1

    mass_of_components = (composition['Hex'] * MASS_HEX +
                          composition['HexNAc'] * MASS_HEXNAC +
                          composition['Neu5Ac'] * MASS_NEU5AC)
    
    mass_lost_from_bonds = num_glycosidic_bonds * MASS_H2O
    
    native_glycan_mass = mass_of_components - mass_lost_from_bonds

    # --- Step 2: Account for amidation ---
    # Reaction: -COOH -> -CONH2. Net change: -O + NH per sialic acid.
    num_amidations = 2
    mass_change_per_amidation = ATOMIC_MASS['N'] + ATOMIC_MASS['H'] - ATOMIC_MASS['O']
    amidation_mass_change = num_amidations * mass_change_per_amidation
    amidated_glycan_mass = native_glycan_mass + amidation_mass_change

    # --- Step 3: Account for permethylation ---
    # Reaction: -H -> -CH3 for each available site. Net change: +CH2 per site.
    # Methylation sites: All -OH groups + N-H of N-acetyl groups.
    # Hexose: 4 -OH; HexNAc: 3 -OH, 1 N-H; Neu5Ac: 4 -OH, 1 N-H
    # The carboxyl groups are now amides and not methylated.
    sites_per_hex = 4
    sites_per_hexnac = 4 # 3 OH + 1 NH
    sites_per_neu5ac = 5 # 4 OH + 1 NH
    
    # Total sites on unlinked monosaccharides
    total_sites_unlinked = (composition['Hex'] * sites_per_hex +
                            composition['HexNAc'] * sites_per_hexnac +
                            composition['Neu5Ac'] * sites_per_neu5ac)
                            
    # Each glycosidic bond consumes two sites (one from each monosaccharide)
    num_methylation_sites = total_sites_unlinked - 2 * num_glycosidic_bonds

    mass_change_per_methylation = ATOMIC_MASS['C'] + 2 * ATOMIC_MASS['H']
    methylation_mass_change = num_methylation_sites * mass_change_per_methylation
    
    final_neutral_mass = amidated_glycan_mass + methylation_mass_change
    
    # --- Step 4: Account for sodiation ---
    sodiated_ion_mass = final_neutral_mass + ATOMIC_MASS['Na']

    # --- Print results ---
    glycan_names = ["A2G(4)2S(3)2", "A2G(4)S(3)S(6)", "A2G(4)2S(6)2"]
    print("The three glycans are isomers and will have the same final m/z.")
    print("\nCalculation steps:")
    print(f"1. Mass of native glycan (Man3GlcNAc4Gal2Sia2): {native_glycan_mass:.4f} Da")
    print(f"2. Mass change from amidation of 2 sialic acids: {amidation_mass_change:.4f} Da")
    print(f"3. Mass of amidated glycan: {amidated_glycan_mass:.4f} Da")
    print(f"4. Number of permethylation sites: {num_methylation_sites}")
    print(f"5. Mass change from permethylation: {methylation_mass_change:.4f} Da")
    print(f"6. Mass of neutral derivatized glycan: {final_neutral_mass:.4f} Da")
    print(f"7. Mass of sodium adduct ([M+Na]+): +{ATOMIC_MASS['Na']:.4f} Da")
    
    print("\n" + "="*40)
    print("Final m/z calculation:")
    print(f"{amidated_glycan_mass:.4f} Da (Amidated Glycan) + {methylation_mass_change:.4f} Da (Permethylation) + {ATOMIC_MASS['Na']:.4f} Da (Sodium) = {sodiated_ion_mass:.4f} Da")
    print("="*40 + "\n")

    for name in glycan_names:
        print(f"The expected m/z for the singly-sodiated ion of {name} is {sodiated_ion_mass:.4f}")

calculate_glycan_mass()

<<<2762.3952>>>