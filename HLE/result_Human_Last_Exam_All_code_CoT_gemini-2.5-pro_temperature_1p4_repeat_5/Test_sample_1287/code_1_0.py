import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a specific derivatized N-glycan.
    The glycan is A2G2S2, which undergoes amidation, permethylation,
    and is detected as a singly sodiated ion.
    """
    # Monoisotopic masses of elements and common adducts
    ATOMIC_MASSES = {
        'H': 1.007825,
        'C': 12.000000,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770,
    }

    # --- Step 1: Calculate the mass of the native glycan ---
    # Monoisotopic masses of monosaccharide residues (loss of H2O)
    mass_H2O = 2 * ATOMIC_MASSES['H'] + ATOMIC_MASSES['O']
    # Hexose (Gal, Man) residue: C6H10O5
    mass_hex_res = 6 * ATOMIC_MASSES['C'] + 10 * ATOMIC_MASSES['H'] + 5 * ATOMIC_MASSES['O']
    # N-acetylhexosamine (GlcNAc) residue: C8H13NO5
    mass_hexnac_res = 8 * ATOMIC_MASSES['C'] + 13 * ATOMIC_MASSES['H'] + ATOMIC_MASSES['N'] + 5 * ATOMIC_MASSES['O']
    # N-acetylneuraminic acid (Neu5Ac) residue: C11H17NO8
    mass_neu5ac_res = 11 * ATOMIC_MASSES['C'] + 17 * ATOMIC_MASSES['H'] + ATOMIC_MASSES['N'] + 8 * ATOMIC_MASSES['O']

    # Glycan composition: 3 Man, 4 GlcNAc, 2 Gal, 2 Neu5Ac
    # Hex = Man + Gal = 5
    # HexNAc = 4
    # Neu5Ac = 2
    native_mass = (5 * mass_hex_res) + (4 * mass_hexnac_res) + (2 * mass_neu5ac_res) + mass_H2O

    # --- Step 2: Calculate mass change from amidation ---
    # The reaction converts -COOH to -CONH2.
    # This is a net change of losing an Oxygen and gaining an NH group.
    amidation_change_per_group = ATOMIC_MASSES['N'] + ATOMIC_MASSES['H'] - ATOMIC_MASSES['O']
    num_sialic_acids = 2
    total_amidation_change = num_sialic_acids * amidation_change_per_group

    # --- Step 3: Calculate mass change from permethylation ---
    # Permethylation replaces H with CH3 on all OH and NH groups.
    # This is a net addition of a CH2 group per site.
    mass_CH2 = ATOMIC_MASSES['C'] + 2 * ATOMIC_MASSES['H']
    
    # Number of methylation sites for an amidated A2G2S2 glycan is 39.
    # This is derived from:
    # - 29 sites on hydroxyl (-OH) groups
    # - 6 sites on N-acetyl (-NH) groups (4 on GlcNAcs, 2 on Neu5Acs)
    # - 4 sites on the two amide (-CONH2) groups (2 per group)
    num_methylation_sites = 39
    total_permethylation_change = num_methylation_sites * mass_CH2

    # --- Step 4: Calculate final m/z for the [M+Na]+ ion ---
    final_mass = native_mass + total_amidation_change + total_permethylation_change
    sodiated_ion_mass = final_mass + ATOMIC_MASSES['Na']

    print("Calculation of the expected m/z for derivatized A2G2S2 glycans:\n")
    print(f"1. Mass of Native Glycan (A2G2S2): {native_mass:.4f} Da")
    print(f"2. Mass Change from Amidation (2x sialic acid): {total_amidation_change:.4f} Da")
    print(f"3. Mass Change from Permethylation (39 sites): +{total_permethylation_change:.4f} Da")
    print(f"4. Mass of Sodium Adduct (Na+): +{ATOMIC_MASSES['Na']:.4f} Da")
    print("-" * 30)

    # Output the final equation with all the numbers
    print("Final m/z Calculation:")
    print(f"m/z = (Native Mass) + (Amidation) + (Permethylation) + (Sodium)")
    print(f"m/z = {native_mass:.4f} + ({total_amidation_change:.4f}) + {total_permethylation_change:.4f} + {ATOMIC_MASSES['Na']:.4f}")
    print(f"m/z = {sodiated_ion_mass:.4f}\n")
    
    print("The three glycans you listed are isomers and do not differ in mass after these modifications.")
    print("Therefore, the expected m/z for the singly-sodiated ion is the same for all three:")
    print(f"Mass for A2G(4)2S(3)2: {sodiated_ion_mass:.4f}")
    print(f"Mass for A2G(4)S(3)S(6): {sodiated_ion_mass:.4f}")
    print(f"Mass for A2G(4)2S(6)2: {sodiated_ion_mass:.4f}")
    
    return sodiated_ion_mass

# Execute the calculation and print the results
calculated_mass = calculate_glycan_mass()
# The final answer in the requested format
print(f"\n<<<{calculated_mass:.4f}>>>")
