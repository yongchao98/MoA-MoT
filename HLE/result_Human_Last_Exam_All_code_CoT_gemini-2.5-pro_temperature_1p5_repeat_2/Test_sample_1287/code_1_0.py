import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized biantennary disialylated N-glycan.
    
    The glycan undergoes two reactions:
    1. Amidation of the two sialic acids' carboxyl groups.
    2. Permethylation of all free hydroxyl and N-acetyl groups.
    
    The final ion is a singly-sodiated adduct [M+Na]+.
    """
    
    # Monoisotopic atomic masses
    atomic_mass = {
        'C': 12.0000000,
        'H': 1.0078250,
        'N': 14.0030740,
        'O': 15.9949150,
        'Na': 22.9897700
    }

    # Elemental composition of individual monosaccharides
    mono_comp = {
        'Neu5Ac':   {'C': 11, 'H': 19, 'N': 1, 'O': 9},
        'Gal':      {'C': 6, 'H': 12, 'N': 0, 'O': 6},
        'GlcNAc':   {'C': 8, 'H': 15, 'N': 1, 'O': 6},
        'Man':      {'C': 6, 'H': 12, 'N': 0, 'O': 6}
    }
    
    water_comp = {'H': 2, 'O': 1}
    
    # Base glycan composition: 2x Neu5Ac, 2x Gal, 4x GlcNAc, 3x Man
    glycan_counts = {'Neu5Ac': 2, 'Gal': 2, 'GlcNAc': 4, 'Man': 3}
    num_monosaccharides = sum(glycan_counts.values())
    num_glycosidic_bonds = num_monosaccharides - 1
    
    # 1. Calculate the initial elemental composition of the glycan
    base_formula = collections.defaultdict(int)
    for mono, count in glycan_counts.items():
        for element, num in mono_comp[mono].items():
            base_formula[element] += num * count
            
    # Subtract water for glycosidic bonds
    for element, num in water_comp.items():
        base_formula[element] -= num * num_glycosidic_bonds
        
    print(f"Step 1: The starting glycan has the elemental formula C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}.")

    # 2. Adjust for amidation of the two sialic acids
    # Reaction: R-COOH -> R-CONH2. Net change: -O +NH
    num_amidations = 2
    amidated_formula = base_formula.copy()
    amidated_formula['O'] -= 1 * num_amidations
    amidated_formula['N'] += 1 * num_amidations
    amidated_formula['H'] += 1 * num_amidations
    
    print(f"Step 2: After amidation of {num_amidations} sialic acids, the formula is C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}.")
    
    # 3. Count methylation sites
    # Sites = all -OH groups + N-H groups of GlcNAc residues.
    # The primary amide N-H from amidation is not methylated under these conditions.
    sites_reducing_glcnac = 4  # C1-OH, C3-OH, C6-OH, N-H
    sites_core_glcnac = 3      # C3-OH, C6-OH, N-H
    sites_branching_man = 2  # C2-OH, C4-OH
    sites_arm_man = 2 + 2      # 2 per arm (3-arm and 6-arm)
    sites_antenna_glcnac = 3 * 2 # 3 per GlcNAc (C3-OH, C6-OH, N-H)
    sites_gal = 3 * 2          # 3 per Galactose residue
    sites_sialic_acid_amide = 5 * 2 # C4,C7,C8,C9-OH + N-H of acetyl
    
    total_methylation_sites = (sites_reducing_glcnac + sites_core_glcnac + sites_branching_man +
                               sites_arm_man + sites_antenna_glcnac + sites_gal + sites_sialic_acid_amide)

    print(f"Step 3: The total number of permethylation sites is {total_methylation_sites}.")

    # 4. Adjust for permethylation
    # Reaction: R-XH -> R-X-CH3. Net change: -H +CH3, which is +CH2
    final_formula = amidated_formula.copy()
    final_formula['C'] += total_methylation_sites
    # H changes by -1 (loss of H) + 3 (gain of CH3) = +2 for each site
    final_formula['H'] += 2 * total_methylation_sites
    
    print(f"Step 4: After permethylation, the final neutral formula is C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}.")

    # 5. Calculate the mass of the neutral molecule
    neutral_mass = 0
    for element, count in final_formula.items():
        neutral_mass += count * atomic_mass[element]
    
    print(f"Step 5: The monoisotopic mass of the neutral derivatized molecule is {neutral_mass:.6f} Da.")

    # 6. Calculate m/z for the singly-sodiated ion [M+Na]+
    sodiated_ion_mass = neutral_mass + atomic_mass['Na']
    # The charge z is +1
    mz_value = sodiated_ion_mass / 1

    print("\n--- Final Result ---")
    print("The three starting glycans are isomers. The amidation and permethylation procedure results in final products that are also isomers of each other.")
    print("Therefore, all three molecules will have the same mass.")
    print(f"Final calculation: Mass({final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}) + Mass(Na) = m/z")
    print(f"{neutral_mass:.6f} + {atomic_mass['Na']:.6f} = {mz_value:.6f}")
    
    print("\nThe expected m/z for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 is:")
    print(f"{mz_value:.4f}")
    
    return mz_value

final_mass = calculate_glycan_mass()
# The final answer format is requested at the very end
# But to show the calculation, I'll print the number. The prompt will take the last line of the script.
# final_answer = f"<<<{final_mass:.4f}>>>"
# print(final_answer) # This would be stripped anyway, so I'll just rely on the final output value.
<<<2734.3525>>>