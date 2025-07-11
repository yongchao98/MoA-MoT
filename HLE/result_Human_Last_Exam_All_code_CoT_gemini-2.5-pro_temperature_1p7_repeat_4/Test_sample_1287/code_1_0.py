import sys
# A function to ensure we don't mix up print statements in the code and the final output
def final_print(text):
    sys.stdout.write(text)

def calculate_glycan_mass():
    """
    Calculates the mass of a permethylated, amidated A2G2S2 glycan ion.
    """
    # Monoisotopic atomic masses from IUPAC
    ATOMIC_MASS = {
        'H': 1.00782503207,
        'C': 12.0000000,
        'N': 14.0030740048,
        'O': 15.99491461956,
        'Na': 22.9897692809,
    }

    # --- Step 1: Calculate the mass of the initial glycan ---
    # Chemical formulas for the monosaccharides
    FORMULA = {
        'HexNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6},  # N-acetylglucosamine (GlcNAc)
        'Hex':    {'C': 6, 'H': 12, 'N': 0, 'O': 6},  # Mannose, Galactose
        'NeuAc':  {'C': 11, 'H': 19, 'N': 1, 'O': 9}, # N-acetylneuraminic acid
        'H2O':    {'C': 0, 'H': 2, 'N': 0, 'O': 1},
    }

    # Composition of the glycan
    COMPOSITION = {
        'HexNAc': 4,
        'Hex': 5, # 3 Mannose + 2 Galactose
        'NeuAc': 2,
    }
    
    # Calculate mass of individual monosaccharides
    mass_hexnac = FORMULA['HexNAc']['C'] * ATOMIC_MASS['C'] + FORMULA['HexNAc']['H'] * ATOMIC_MASS['H'] + FORMULA['HexNAc']['N'] * ATOMIC_MASS['N'] + FORMULA['HexNAc']['O'] * ATOMIC_MASS['O']
    mass_hex = FORMULA['Hex']['C'] * ATOMIC_MASS['C'] + FORMULA['Hex']['H'] * ATOMIC_MASS['H'] + FORMULA['Hex']['N'] * ATOMIC_MASS['N'] + FORMULA['Hex']['O'] * ATOMIC_MASS['O']
    mass_neuac = FORMULA['NeuAc']['C'] * ATOMIC_MASS['C'] + FORMULA['NeuAc']['H'] * ATOMIC_MASS['H'] + FORMULA['NeuAc']['N'] * ATOMIC_MASS['N'] + FORMULA['NeuAc']['O'] * ATOMIC_MASS['O']
    mass_h2o = FORMULA['H2O']['H'] * ATOMIC_MASS['H'] + FORMULA['H2O']['O'] * ATOMIC_MASS['O']

    # Total number of monosaccharides is 4+5+2 = 11. This requires 10 glycosidic bonds.
    num_linkages = sum(COMPOSITION.values()) - 1
    
    mass_of_all_sugars = (COMPOSITION['HexNAc'] * mass_hexnac) + (COMPOSITION['Hex'] * mass_hex) + (COMPOSITION['NeuAc'] * mass_neuac)
    initial_glycan_mass = mass_of_all_sugars - (num_linkages * mass_h2o)

    # --- Step 2: Calculate mass change from amidation ---
    # The reaction is R-COOH -> R-CONH2. The net change is the removal of one O and the addition of one N and one H.
    amidation_mass_change_per_site = ATOMIC_MASS['N'] + ATOMIC_MASS['H'] - ATOMIC_MASS['O']
    num_amidation_sites = 2 # One for each sialic acid
    total_amidation_change = num_amidation_sites * amidation_mass_change_per_site
    amidated_glycan_mass = initial_glycan_mass + total_amidation_change

    # --- Step 3: Calculate mass change from permethylation ---
    # Permethylation replaces an H with a CH3 group, a net addition of CH2.
    permethylation_mass_change_per_site = ATOMIC_MASS['C'] + 2 * ATOMIC_MASS['H']
    
    # Count methylation sites on the amidated glycan:
    # 34 sites on the base glycan skeleton (from all -OH and N-acetyl -NH groups).
    # 2 sites per -CONH2 group formed, and there are 2 such groups. (2 * 2 = 4).
    num_permethylation_sites = 34 + 4
    total_permethylation_change = num_permethylation_sites * permethylation_mass_change_per_site

    # --- Step 4: Calculate final ion mass ---
    final_neutral_mass = amidated_glycan_mass + total_permethylation_change
    sodiated_ion_mass = final_neutral_mass + ATOMIC_MASS['Na']

    # --- Step 5: Print the explanation and result ---
    final_print("Since A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers with the same number of reactive sites, they will all have the same final mass.\n\n")
    final_print("The calculation for the final m/z is as follows:\n\n")
    final_print("Final Mass = (Mass of Initial Glycan) + (Mass Change from Amidation) + (Mass Change from Permethylation) + (Mass of Sodium)\n\n")
    
    final_print(f"Mass of Initial Glycan:\n")
    final_print(f"  ({COMPOSITION['HexNAc']} * Mass(HexNAc)) + ({COMPOSITION['Hex']} * Mass(Hex)) + ({COMPOSITION['NeuAc']} * Mass(NeuAc)) - ({num_linkages} * Mass(H2O))\n")
    final_print(f"= ({COMPOSITION['HexNAc']} * {mass_hexnac:.4f}) + ({COMPOSITION['Hex']} * {mass_hex:.4f}) + ({COMPOSITION['NeuAc']} * {mass_neuac:.4f}) - ({num_linkages} * {mass_h2o:.4f})\n")
    final_print(f"= {initial_glycan_mass:.4f} Da\n\n")

    final_print(f"Mass Change from Amidation ({num_amidation_sites} sites):\n")
    final_print(f"  {num_amidation_sites} * (Mass(N) + Mass(H) - Mass(O))\n")
    final_print(f"= {num_amidation_sites} * ({amidation_mass_change_per_site:.4f})\n")
    final_print(f"= {total_amidation_change:.4f} Da\n\n")

    final_print(f"Mass Change from Permethylation ({num_permethylation_sites} sites):\n")
    final_print(f"  {num_permethylation_sites} * (Mass(C) + 2 * Mass(H))\n")
    final_print(f"= {num_permethylation_sites} * ({permethylation_mass_change_per_site:.4f})\n")
    final_print(f"= {total_permethylation_change:.4f} Da\n\n")

    final_print(f"Mass of Sodium Ion:\n= {ATOMIC_MASS['Na']:.4f} Da\n\n")

    final_print(f"Final Observed Mass [M+Na]+:\n")
    final_print(f"= {initial_glycan_mass:.4f} + {total_amidation_change:.4f} + {total_permethylation_change:.4f} + {ATOMIC_MASS['Na']:.4f}\n")
    final_print(f"= {sodiated_ion_mass:.4f} m/z\n")
    
    return sodiated_ion_mass

observed_mass = calculate_glycan_mass()
final_print(f"\nTherefore, the mass observed for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 should all be {observed_mass:.4f} m/z.")
# Use a special format for the final numerical answer.
final_print(f"\n<<<{observed_mass:.4f}>>>")