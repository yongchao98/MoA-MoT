def calculate_glycan_mass():
    """
    Calculates the expected m/z for permethylated and amidated A2G2S2 glycans.
    """
    # Monoisotopic atomic masses
    ATOMIC_MASS = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989769
    }

    # Elemental composition of monosaccharides
    MONOSACCHARIDE_FORMULA = {
        'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6},
        'Hexose': {'C': 6, 'H': 12, 'O': 6},  # Man and Gal are hexoses
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9}
    }
    
    # Water formula for calculating glycosidic bond formation
    H2O_FORMULA = {'H': 2, 'O': 1}

    # Composition of the A2G2S2 glycan
    glycan_composition = {
        'GlcNAc': 2,
        'Hexose': 5, # 3 Mannose + 2 Galactose
        'Neu5Ac': 2
    }
    
    glycan_names = ["A2G(4)2S(3)2", "A2G(4)S(3)S(6)", "A2G(4)2S(6)2"]

    # --- Step 1: Calculate the formula of the native glycan ---
    base_formula = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    num_monosaccharides = 0
    for mono, count in glycan_composition.items():
        num_monosaccharides += count
        for atom, atom_count in MONOSACCHARIDE_FORMULA[mono].items():
            base_formula[atom] += atom_count * count

    # Subtract water molecules for glycosidic bonds (n-1 bonds for n units)
    num_bonds = num_monosaccharides - 1
    base_formula['H'] -= H2O_FORMULA['H'] * num_bonds
    base_formula['O'] -= H2O_FORMULA['O'] * num_bonds

    # --- Step 2: Adjust formula for amidation ---
    # Reaction on two sialic acids: R-COOH -> R-CONH2. Net change per reaction: -O +NH
    amidated_formula = base_formula.copy()
    num_sialic_acids = glycan_composition['Neu5Ac']
    amidated_formula['O'] -= 1 * num_sialic_acids
    amidated_formula['N'] += 1 * num_sialic_acids
    amidated_formula['H'] += 1 * num_sialic_acids

    # --- Step 3: Adjust formula for permethylation ---
    # Count methylation sites (all -OH and -NH groups)
    # GlcNAc: 3 -OH + 1 -NH = 4 sites
    # Hexose: 4 -OH = 4 sites
    # Neu5Ac: 3 -OH + 1 -NH = 4 sites (COOH is now CONH2, which is not methylated)
    # Extra site for free reducing end: 1 site
    sites_per_mono = {
        'GlcNAc': 4,
        'Hexose': 4,
        'Neu5Ac': 4
    }
    total_methylation_sites = 0
    for mono, count in glycan_composition.items():
        total_methylation_sites += sites_per_mono[mono] * count
    total_methylation_sites += 1 # For the free reducing end

    # Add CH2 for each methylation site (replaces H with CH3)
    final_formula = amidated_formula.copy()
    final_formula['C'] += 1 * total_methylation_sites
    final_formula['H'] += 2 * total_methylation_sites
    
    # --- Step 4: Calculate the final mass of the [M+Na]+ ion ---
    mass_M = (final_formula['C'] * ATOMIC_MASS['C'] +
              final_formula['H'] * ATOMIC_MASS['H'] +
              final_formula['N'] * ATOMIC_MASS['N'] +
              final_formula['O'] * ATOMIC_MASS['O'])
    
    mass_M_plus_Na = mass_M + ATOMIC_MASS['Na']

    # --- Step 5: Print the results ---
    print("The three glycans are isomers and will have the same mass after the described chemical modifications.")
    print("-" * 30)
    
    print(f"The final elemental formula of the modified glycan (M) is: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")
    print(f"The observed ion is [M+Na]+.")
    print("\nFinal Mass Calculation:")
    print(f"Mass = (C atoms * C_mass) + (H atoms * H_mass) + (N atoms * N_mass) + (O atoms * O_mass) + Na_mass")
    print(f"Mass = ({final_formula['C']} * {ATOMIC_MASS['C']}) + ({final_formula['H']} * {ATOMIC_MASS['H']}) + ({final_formula['N']} * {ATOMIC_MASS['N']}) + ({final_formula['O']} * {ATOMIC_MASS['O']}) + {ATOMIC_MASS['Na']}")
    
    print("\nExpected m/z for [M+Na]+:")
    for name in glycan_names:
        print(f"{name}: {mass_M_plus_Na:.4f}")

calculate_glycan_mass()