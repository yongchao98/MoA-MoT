def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized (amidated and permethylated)
    singly-sodiated A2G2S2 glycan.
    """
    # Monoisotopic masses of the elements and sodium
    ATOMIC_MASSES = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Step 1: Determine the molecular formula of the native glycan (A2G2S2)
    # Monosaccharide compositions
    glycan_composition = {
        'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6, 'count': 4},
        'Man':    {'C': 6, 'H': 12, 'N': 0, 'O': 6, 'count': 3},
        'Gal':    {'C': 6, 'H': 12, 'N': 0, 'O': 6, 'count': 2},
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9, 'count': 2},
    }
    
    # Sum of all atoms before forming the glycan
    total_atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    num_residues = 0
    for residue, data in glycan_composition.items():
        num_residues += data['count']
        for atom, count in data.items():
            if atom in total_atoms:
                total_atoms[atom] += count * data['count']

    # Subtract water molecules (H2O) for each glycosidic bond (n-1 bonds for n residues)
    num_bonds = num_residues - 1
    native_formula = total_atoms.copy()
    native_formula['H'] -= 2 * num_bonds
    native_formula['O'] -= 1 * num_bonds

    # Step 2: Apply amidation reaction
    # Replaces -COOH with -CONH2 for each of the 2 sialic acids.
    # Net change per reaction: replace OH with NH2.
    # Formula change: -O, -H, +N, +2H  =>  -O, +H, +N
    num_amidations = 2
    amidated_formula = native_formula.copy()
    amidated_formula['O'] -= num_amidations
    amidated_formula['H'] += num_amidations
    amidated_formula['N'] += num_amidations

    # Step 3: Determine the number of permethylation sites (P)
    # Sum of all available -OH and -NH groups on the amidated glycan
    sites = {
        "Reducing end GlcNAc (linked at C4)": 4, # 3 OH + 1 NH
        "Core GlcNAc (linked at C1, C4)": 3, # 2 OH + 1 NH
        "Core Man (linked at C1, C3, C6)": 0, # 0 OH
        "Branch point a1-3 Man (linked at C1, C2)": 2, # 2 OH
        "Branch point a1-6 Man (linked at C1, C2)": 2, # 2 OH
        "Antenna GlcNAc x2 (linked at C1, C4)": 2 * 3, # 2 * (2 OH + 1 NH)
        "Antenna Gal x2 (linked at C1, C3/C6)": 2 * 3, # 2 * (3 OH)
        "Terminal Amido-Neu5Ac x2 (linked at C2)": 2 * 7 # 2 * (4 OH + 1 N-H_acetyl + 2 N-H_amide)
    }
    num_methylation_sites = sum(sites.values())
    
    # Step 4: Apply permethylation reaction
    # Replaces H with CH3 for each site.
    # Net change per reaction: +CH2
    permethylated_formula = amidated_formula.copy()
    permethylated_formula['C'] += num_methylation_sites
    permethylated_formula['H'] += 2 * num_methylation_sites

    # Step 5: Calculate the final m/z for the [M+Na]+ ion
    # Mass of the neutral derivatized molecule
    C = permethylated_formula['C']
    H = permethylated_formula['H']
    N = permethylated_formula['N']
    O = permethylated_formula['O']

    mass_C = C * ATOMIC_MASSES['C']
    mass_H = H * ATOMIC_MASSES['H']
    mass_N = N * ATOMIC_MASSES['N']
    mass_O = O * ATOMIC_MASSES['O']
    
    neutral_mass = mass_C + mass_H + mass_N + mass_O
    
    # Add mass of sodium for the sodiated ion
    sodiated_mz = neutral_mass + ATOMIC_MASSES['Na']

    print("The three glycans are isomers. After the specified derivatization, they all have the same molecular formula and therefore the same mass.")
    print("The calculation for the expected m/z is as follows:\n")
    print(f"1. Final Molecular Formula: C{C} H{H} N{N} O{O}")
    print(f"2. Number of Methylation Sites: {num_methylation_sites}")
    print(f"3. Mass Calculation for the [M+Na]+ ion:")
    print(f"   m/z = Mass(C{C} H{H} N{N} O{O}) + Mass(Na)")
    print(f"   m/z = ({C} * {ATOMIC_MASSES['C']:.6f}) + ({H} * {ATOMIC_MASSES['H']:.6f}) + ({N} * {ATOMIC_MASSES['N']:.6f}) + ({O} * {ATOMIC_MASSES['O']:.6f}) + {ATOMIC_MASSES['Na']:.6f}")
    print(f"   m/z = {mass_C:.6f} + {mass_H:.6f} + {mass_N:.6f} + {mass_O:.6f} + {ATOMIC_MASSES['Na']:.6f}")
    print(f"   m/z = {neutral_mass:.6f} + {ATOMIC_MASSES['Na']:.6f}")
    print(f"   m/z = {sodiated_mz:.4f}\n")
    print("Therefore, the same mass should be observed for all three glycans.")

calculate_glycan_mass()
print(f'<<<2762.3838>>>')