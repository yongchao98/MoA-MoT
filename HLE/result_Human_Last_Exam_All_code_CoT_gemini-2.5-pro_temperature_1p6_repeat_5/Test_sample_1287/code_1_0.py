def calculate_glycan_mass():
    """
    Calculates the m/z of a permethylated, amidated, sodiated A2G2S2 glycan.
    """
    # Monoisotopic masses of the most common isotopes
    atom_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
        'e-': 0.000549
    }

    # Step 1: Define the initial formula of the A2G2S2 glycan
    # Composition: 5 Hexose (Man/Gal), 4 GlcNAc, 2 Neu5Ac
    # Sum of residues: 5*C6H10O5 + 4*C8H13NO5 + 2*C11H17NO8 + H2O (for reducing end)
    # Total Formula: C84 H138 N6 O62
    native_glycan_formula = {'C': 84, 'H': 138, 'N': 6, 'O': 62}

    # Step 2: Account for the amidation of two sialic acids
    # Reaction: R-COOH -> R-CONH2. This is a net change of -OH + NH2 per reaction.
    # Change in formula for two reactions: -O2H2 + N2H4 = N2 H2 O-2
    amidated_glycan_formula = native_glycan_formula.copy()
    amidated_glycan_formula['O'] -= 2
    amidated_glycan_formula['N'] += 2
    amidated_glycan_formula['H'] += 2 # -2H from OH, +4H from NH2

    # Step 3: Account for permethylation
    # We methylate all -OH and -NH protons.
    # Number of sites:
    # - OH groups in native A2G2S2: 21
    # - NH groups in native A2G2S2 (from N-acetyl groups): 6 (4 from GlcNAc, 2 from Neu5Ac)
    # - New NH groups from amidation (2 amides, each -CONH2): 2 * 2 = 4
    # Total methylation sites = 21 + 6 + 4 = 31
    num_methylations = 31

    # For each methylation, an H is replaced by a CH3 group.
    # Net change per methylation: +C +H2
    final_neutral_formula = amidated_glycan_formula.copy()
    final_neutral_formula['C'] += num_methylations
    final_neutral_formula['H'] += num_methylations * 2

    # Step 4: Calculate the monoisotopic mass of the final neutral molecule
    m_neutral = (final_neutral_formula['C'] * atom_mass['C'] +
                 final_neutral_formula['H'] * atom_mass['H'] +
                 final_neutral_formula['N'] * atom_mass['N'] +
                 final_neutral_formula['O'] * atom_mass['O'])

    # Step 5: Calculate the m/z for the singly-sodiated ion [M+Na]+
    # Mass of sodium ion = Mass of sodium atom - mass of electron
    m_na_ion = atom_mass['Na'] - atom_mass['e-']
    final_mz = m_neutral + m_na_ion

    # --- Output Results ---
    print("The three glycans are isomers and will have the identical mass after the specified reactions.\n")
    print("The calculation steps are as follows:")
    print(f"1. Final chemical formula of the modified neutral glycan: C{final_neutral_formula['C']}H{final_neutral_formula['H']}N{final_neutral_formula['N']}O{final_neutral_formula['O']}")
    
    print("\n2. Mass calculation for the sodiated ion [M+Na]+:")
    
    # Print the full equation for transparency
    print(f"m/z = (C * {atom_mass['C']}) + (H * {atom_mass['H']}) + (N * {atom_mass['N']}) + (O * {atom_mass['O']}) + (Na+)")
    print(f"m/z = ({final_neutral_formula['C']} * {atom_mass['C']}) + ({final_neutral_formula['H']} * {atom_mass['H']}) + "
          f"({final_neutral_formula['N']} * {atom_mass['N']}) + ({final_neutral_formula['O']} * {atom_mass['O']}) + "
          f"{m_na_ion:.6f}")
    
    c_mass_contrib = final_neutral_formula['C'] * atom_mass['C']
    h_mass_contrib = final_neutral_formula['H'] * atom_mass['H']
    n_mass_contrib = final_neutral_formula['N'] * atom_mass['N']
    o_mass_contrib = final_neutral_formula['O'] * atom_mass['O']

    print(f"m/z = {c_mass_contrib:.4f} + {h_mass_contrib:.4f} + {n_mass_contrib:.4f} + {o_mass_contrib:.4f} + {m_na_ion:.4f}")
    print(f"m/z = {final_mz:.4f}")

    print("\nTherefore, the expected mass for all three of your glycan preparations is:")
    print(f"A2G(4)2S(3)2: m/z = {final_mz:.4f}")
    print(f"A2G(4)S(3)S(6): m/z = {final_mz:.4f}")
    print(f"A2G(4)2S(6)2: m/z = {final_mz:.4f}")

calculate_glycan_mass()