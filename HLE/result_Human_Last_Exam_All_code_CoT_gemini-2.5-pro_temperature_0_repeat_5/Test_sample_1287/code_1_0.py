import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a specific modified and permethylated glycan.
    """
    # Monoisotopic atomic masses
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770
    }

    # Step 1: Define the base composition and formula of the A2G2S2 glycan
    # Composition: 4 GlcNAc, 3 Man, 2 Gal, 2 Neu5Ac
    monosaccharide_formulas = {
        'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6},
        'Man':    {'C': 6, 'H': 12, 'O': 6},
        'Gal':    {'C': 6, 'H': 12, 'O': 6},
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9}
    }
    glycan_composition = {
        'GlcNAc': 4,
        'Man': 3,
        'Gal': 2,
        'Neu5Ac': 2
    }

    # Sum the atoms from all monosaccharides
    total_atoms = collections.defaultdict(int)
    for saccharide, count in glycan_composition.items():
        for atom, num in monosaccharide_formulas[saccharide].items():
            total_atoms[atom] += count * num

    # Subtract water molecules for glycosidic bonds
    # Number of residues = 4 + 3 + 2 + 2 = 11
    # Number of bonds = 11 - 1 = 10
    num_bonds = sum(glycan_composition.values()) - 1
    total_atoms['H'] -= 2 * num_bonds
    total_atoms['O'] -= 1 * num_bonds
    
    base_formula = dict(total_atoms)

    # Step 2: Apply amidation of the two sialic acids
    # Reaction: -COOH -> -CONH2, which is a net change of -O +NH per site
    amidation_sites = 2
    total_atoms['O'] -= 1 * amidation_sites
    total_atoms['N'] += 1 * amidation_sites
    total_atoms['H'] += 1 * amidation_sites
    
    amidated_formula = dict(total_atoms)

    # Step 3: Apply exhaustive permethylation
    # This adds a methyl group (CH3) for every H on an O or N atom, a net change of +CH2
    # Count all free -OH and -NH sites on the amidated glycan
    # Free -OH sites: 31 (30 on the glycan body + 1 at the reducing end)
    # Free N-acetyl -NH sites: 6 (4 from GlcNAc, 2 from Neu5Ac)
    # Free primary amide -NH2 sites: 2 (from the two amidated Neu5Ac, each has 2 N-H bonds)
    oh_sites = 31
    n_acetyl_nh_sites = 6
    primary_amide_nh_sites = 2 * 2 # 2 amides, 2 N-H each
    total_methylation_sites = oh_sites + n_acetyl_nh_sites + primary_amide_nh_sites
    
    # Add atoms for methylation (+CH2 per site)
    total_atoms['C'] += 1 * total_methylation_sites
    total_atoms['H'] += 2 * total_methylation_sites
    
    final_formula = dict(total_atoms)

    # Step 4 & 5: Calculate the mass of the final molecule [M+Na]+
    mass_m = (final_formula['C'] * atomic_mass['C'] +
              final_formula['H'] * atomic_mass['H'] +
              final_formula['N'] * atomic_mass['N'] +
              final_formula['O'] * atomic_mass['O'])
              
    mass_m_na = mass_m + atomic_mass['Na']

    # --- Output the results ---
    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("After amidation and permethylation, they will all have the same mass.\n")
    print(f"The final chemical formula of the modified glycan is C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}.\n")
    print("The mass of the singly-sodiated ion [M+Na]+ is calculated as follows:")
    
    c_mass_str = f"({final_formula['C']} * {atomic_mass['C']:.6f})"
    h_mass_str = f"({final_formula['H']} * {atomic_mass['H']:.6f})"
    n_mass_str = f"({final_formula['N']} * {atomic_mass['N']:.6f})"
    o_mass_str = f"({final_formula['O']} * {atomic_mass['O']:.6f})"
    na_mass_str = f"{atomic_mass['Na']:.6f}"
    
    print(f"Mass = C_mass + H_mass + N_mass + O_mass + Na_mass")
    print(f"Mass = {c_mass_str} + {h_mass_str} + {n_mass_str} + {o_mass_str} + {na_mass_str}\n")
    
    c_val = final_formula['C'] * atomic_mass['C']
    h_val = final_formula['H'] * atomic_mass['H']
    n_val = final_formula['N'] * atomic_mass['N']
    o_val = final_formula['O'] * atomic_mass['O']
    
    print(f"Mass = {c_val:.6f} + {h_val:.6f} + {n_val:.6f} + {o_val:.6f} + {atomic_mass['Na']:.6f}\n")
    
    print(f"The expected m/z for the singly-sodiated ions is {mass_m_na:.4f}.")

calculate_glycan_mass()
<<<2818.4466>>>