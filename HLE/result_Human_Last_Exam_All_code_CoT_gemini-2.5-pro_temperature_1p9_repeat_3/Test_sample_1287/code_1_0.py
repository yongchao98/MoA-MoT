def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized biantennary glycan.
    The glycan is A2G2S2, which has been amidated and then permethylated.
    The final ion is a singly-sodiated species [M+Na]+.
    """

    # Monoisotopic masses of elements and sodium
    atomic_mass = {
        'H': 1.007825,
        'C': 12.000000,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770
    }

    # Formulas of individual free monosaccharides
    mono_formulas = {
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9},
        'Gal':    {'C': 6, 'H': 12, 'N': 0, 'O': 6},
        'Man':    {'C': 6, 'H': 12, 'N': 0, 'O': 6},
        'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6}
    }

    # Composition of the glycan
    composition = {
        'Neu5Ac': 2,
        'Gal': 2,
        'Man': 3,
        'GlcNAc': 4
    }
    
    total_monos = sum(composition.values())
    num_links = total_monos - 1

    # 1. Calculate the formula of the intact native glycan
    native_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for mono, count in composition.items():
        for element, num in mono_formulas[mono].items():
            native_formula[element] += num * count
    
    # Subtract H2O for each glycosidic linkage
    native_formula['H'] -= 2 * num_links
    native_formula['O'] -= 1 * num_links

    # 2. Adjust formula for amidation of the two sialic acids
    # Reaction: R-COOH -> R-CONH2. Per reaction: -O, +N, +H
    amidated_formula = native_formula.copy()
    num_sialic_acids = composition['Neu5Ac']
    amidated_formula['O'] -= 1 * num_sialic_acids
    amidated_formula['N'] += 1 * num_sialic_acids
    amidated_formula['H'] += 1 * num_sialic_acids

    # 3. Calculate number of methylation sites and adjust formula
    # Sites = (all OH groups) + (all GlcNAc NH groups) + (all new amide NH2 groups)
    
    # Count OH groups on free monosaccharides
    # Neu5Ac: 5, Gal: 5, Man: 5, GlcNAc: 4 (incl. anomeric)
    total_oh_free = (composition['Neu5Ac'] * 5 + 
                     composition['Gal'] * 5 + 
                     composition['Man'] * 5 + 
                     composition['GlcNAc'] * 4)

    # Free OHs in glycan = Total OHs - 2 per linkage
    glycan_oh_sites = total_oh_free - 2 * num_links

    # NH sites from GlcNAc
    glcnac_nh_sites = composition['GlcNAc']

    # NH2 sites from amidated sialic acids (2 N-H bonds per amide)
    amide_nh_sites = num_sialic_acids * 2

    total_methylation_sites = glycan_oh_sites + glcnac_nh_sites + amide_nh_sites

    # Adjust formula for permethylation
    # Per site: -H, +CH3 (net gain of CH2)
    final_formula = amidated_formula.copy()
    final_formula['C'] += total_methylation_sites
    final_formula['H'] += 2 * total_methylation_sites

    # 4. Calculate the final mass
    mass_M = (final_formula['C'] * atomic_mass['C'] +
              final_formula['H'] * atomic_mass['H'] +
              final_formula['N'] * atomic_mass['N'] +
              final_formula['O'] * atomic_mass['O'])
    
    mass_M_plus_Na = mass_M + atomic_mass['Na']

    # --- Output the results ---
    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("After amidation and permethylation, they remain isomers and will have the same mass.")
    print("\n--- Calculation Steps ---")
    print(f"1. Base Glycan Composition: {composition}")
    print(f"2. Chemical Reactions: Amidation of 2 Sialic Acids, then Permethylation.")
    print(f"3. Number of Methylation Sites: {total_methylation_sites}")
    print(f"4. Final Molecular Formula: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")
    print(f"5. Ion Detected: [M+Na]+")
    print("\n--- Final Mass Calculation ---")
    print("Mass Equation:")
    print(f"m/z = (C * {atomic_mass['C']}) + (H * {atomic_mass['H']}) + (N * {atomic_mass['N']}) + (O * {atomic_mass['O']}) + Na")
    
    print(f"\nm/z = ({final_formula['C']} * {atomic_mass['C']}) + ({final_formula['H']} * {atomic_mass['H']:.6f}) + "
          f"({final_formula['N']} * {atomic_mass['N']:.6f}) + ({final_formula['O']} * {atomic_mass['O']:.6f}) + {atomic_mass['Na']:.6f}")
          
    print(f"\nm/z = {final_formula['C'] * atomic_mass['C']:.4f} + {final_formula['H'] * atomic_mass['H']:.4f} + "
          f"{final_formula['N'] * atomic_mass['N']:.4f} + {final_formula['O'] * atomic_mass['O']:.4f} + {atomic_mass['Na']:.4f}")

    print(f"\nCalculated neutral mass M = {mass_M:.4f} Da")
    print(f"Calculated m/z for [M+Na]+ = {mass_M_plus_Na:.4f}")
    
    print("\nTherefore, the expected mass for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 is:")
    print(f"{mass_M_plus_Na:.4f}")


calculate_glycan_mass()
