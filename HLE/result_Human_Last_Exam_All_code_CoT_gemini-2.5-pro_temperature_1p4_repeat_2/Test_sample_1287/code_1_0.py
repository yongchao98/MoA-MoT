import collections

def calculate_glycan_mz():
    """
    Calculates the m/z of a permethylated and amidated A2G2S2 glycan.
    """
    # Monoisotopic atomic masses
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Elemental composition of monosaccharide residues (after forming glycosidic bond)
    # Hexose includes Man and Gal
    residue_composition = {
        'Hex': {'C': 6, 'H': 10, 'O': 5},
        'GlcNAc': {'C': 8, 'H': 13, 'N': 1, 'O': 5},
        'Neu5Ac': {'C': 11, 'H': 17, 'N': 1, 'O': 8},
    }

    # Glycan structure: Man(3)GlcNAc(4)Gal(2)Neu5Ac(2)
    # A biantennary glycan with two galactose and two sialic acids.
    # We treat Galactose as a Hexose for composition counting.
    glycan_counts = {'Hex': 5, 'GlcNAc': 4, 'Neu5Ac': 2} # 3 Man + 2 Gal = 5 Hex

    # --- Step 1: Calculate the base formula of the free glycan ---
    base_formula = collections.defaultdict(int)
    for residue, count in glycan_counts.items():
        for atom, num in residue_composition[residue].items():
            base_formula[atom] += num * count

    # Add back one water molecule for the free reducing end
    base_formula['H'] += 2
    base_formula['O'] += 1

    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("They will all have the same mass after the described derivatization.\n")
    print("--- Calculation Steps ---")
    print(f"1. Base glycan chemical formula: C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}")

    # --- Step 2: Model the amidation reaction ---
    # Each of the 2 sialic acids has a -COOH group converted to -CONH2.
    # This replaces an -OH group with an -NH2 group.
    # Net change per reaction: -O, -H, +N, +2H  =>  +N, +H, -O
    # Total change for 2 sialic acids: +2N, +2H, -2O
    amidated_formula = base_formula.copy()
    amidated_formula['N'] += 2
    amidated_formula['H'] += 2
    amidated_formula['O'] -= 2
    print(f"2. Formula after amidation of 2 sialic acids: C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}")


    # --- Step 3: Model the permethylation reaction ---
    # We count all acidic protons (-OH and -NH) that get replaced by methyl groups.
    # A) -OH groups:
    #    - Reducing GlcNAc: 4 OH (including anomeric)
    #    - Core GlcNAc: 2 OH
    #    - Core Man (branch): 2 OH
    #    - Arm Man x2: 2 * 2 = 4 OH
    #    - Outer GlcNAc x2: 2 * 2 = 4 OH
    #    - Gal x2: 2 * 3 = 6 OH
    #    - Neu5Ac x2: 2 * 4 = 8 OH
    #    Total OH groups = 4+2+2+4+4+6+8 = 30
    num_oh_sites = 30

    # B) -NH groups:
    #    - GlcNAc N-acetyl groups: 4
    #    - Neu5Ac N-acetyl groups: 2
    #    - Primary amides from amidation (-CONH2): 2 amides * 2 protons/amide = 4
    #    Total NH groups = 4+2+4 = 10
    num_nh_sites = 10

    total_methylation_sites = num_oh_sites + num_nh_sites
    print(f"3. Number of sites for permethylation (-OH and -NH): {total_methylation_sites}")

    # For each site, one -H is replaced by a -CH3 group. Net change is +CH2.
    # Total change: +(total_sites) C, +(total_sites * 2) H
    final_formula = amidated_formula.copy()
    final_formula['C'] += total_methylation_sites
    final_formula['H'] += total_methylation_sites * 2
    print(f"4. Final formula after permethylation: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")


    # --- Step 4 & 5: Calculate final mass and m/z ---
    mass_final_glycan = (final_formula['C'] * atomic_mass['C'] +
                         final_formula['H'] * atomic_mass['H'] +
                         final_formula['N'] * atomic_mass['N'] +
                         final_formula['O'] * atomic_mass['O'])

    mass_na = atomic_mass['Na']
    charge = 1
    final_mz = (mass_final_glycan + mass_na) / charge
    
    print("\n--- Final Result ---")
    print("The calculated m/z is the same for all three glycans.")
    # Printing the final equation with all its components as requested
    print(f"\nFinal Equation: m/z = (Mass_Derivatized_Glycan + Mass_Na) / Charge")
    print(f"m/z = ({mass_final_glycan:.4f} + {mass_na:.4f}) / {charge}")
    print(f"\nExpected m/z for [M+Na]+ ion:")
    print(f"{final_mz:.4f}")

    return final_mz

# Execute the calculation and print the final answer
final_answer = calculate_glycan_mz()
# The final answer in the required format
print(f"\n<<<{final_answer:.4f}>>>")