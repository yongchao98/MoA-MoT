import collections

def calculate_mass(formula_dict, atomic_masses_dict):
    """Calculates the monoisotopic mass from a chemical formula dictionary."""
    mass = 0.0
    for atom, count in formula_dict.items():
        mass += count * atomic_masses_dict[atom]
    return mass

def main():
    """
    Calculates the m/z of derivatized A2G2S2 glycans.
    """
    # 1. Define Monoisotopic Atomic Masses
    ATOMIC_MASSES = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # 2. Define Chemical Formulas for Residues and Molecules
    RESIDUE_FORMULAS = {
        'Hex':    {'C': 6, 'H': 10, 'O': 5}, # Galactose, Mannose
        'GlcNAc': {'C': 8, 'H': 13, 'N': 1, 'O': 5},
        'Neu5Ac': {'C': 11, 'H': 17, 'N': 1, 'O': 8},
        'H2O':    {'C': 0, 'H': 2, 'O': 1},
    }

    # Calculate masses of residues
    mass_hex_res = calculate_mass(RESIDUE_FORMULAS['Hex'], ATOMIC_MASSES)
    mass_glcnac_res = calculate_mass(RESIDUE_FORMULAS['GlcNAc'], ATOMIC_MASSES)
    mass_neu5ac_res = calculate_mass(RESIDUE_FORMULAS['Neu5Ac'], ATOMIC_MASSES)
    mass_h2o = calculate_mass(RESIDUE_FORMULAS['H2O'], ATOMIC_MASSES)

    # 3. Calculate the mass of the native A2G2S2 glycan
    # Composition: 5 Hexose (3 Man, 2 Gal), 4 GlcNAc, 2 Neu5Ac
    # Formula = 5*Hex_res + 4*GlcNAc_res + 2*Neu5Ac_res + 1*H2O (for reducing end)
    mass_native_glycan = (5 * mass_hex_res) + \
                         (4 * mass_glcnac_res) + \
                         (2 * mass_neu5ac_res) + \
                         mass_h2o

    # 4. Calculate mass change from amidation of two sialic acids
    # Reaction: -COOH -> -CONH2. Net change: -O +NH
    mass_change_amidation = 2 * (ATOMIC_MASSES['N'] + ATOMIC_MASSES['H'] - ATOMIC_MASSES['O'])
    
    # 5. Calculate mass change from permethylation
    # Number of methylation sites for N-glycan = 3*H + 3*N + 5*S - 3
    # H=Hexose=5, N=HexNAc=4, S=Sialic Acid=2
    num_sites = (3 * 5) + (3 * 4) + (5 * 2) - 3  # 15 + 12 + 10 - 3 = 34
    
    # Mass change per site: H replaced by CH3. Net change: +CH2
    mass_change_per_methylation = ATOMIC_MASSES['C'] + 2 * ATOMIC_MASSES['H']
    total_mass_change_methylation = num_sites * mass_change_per_methylation

    # 6. Calculate mass of the final derivatized glycan (M)
    mass_final_glycan = mass_native_glycan + mass_change_amidation + total_mass_change_methylation

    # 7. Calculate m/z of the singly sodiated ion [M+Na]+
    mass_na = ATOMIC_MASSES['Na']
    final_mz = mass_final_glycan + mass_na

    # Print the results and explanation
    print("Calculation of the expected m/z for derivatized sialylated glycans")
    print("-" * 70)
    print("The three specified glycans (A2G(4)2S(3)2, A2G(4)S(3)S(6), A2G(4)2S(6)2) are isomers.")
    print("The derivatization reactions (amidation and permethylation) are not linkage-specific.")
    print("Therefore, all three will have the same theoretical mass-to-charge ratio (m/z).\n")

    print("Step-by-step mass calculation:")
    print(f"1. Mass of native A2G2S2 glycan:                             {mass_native_glycan:10.4f} Da")
    print(f"2. Mass change from amidation of 2 sialic acids (-O +NH)*2:  {mass_change_amidation:10.4f} Da")
    print(f"3. Mass change from permethylation ({num_sites} sites * CH2):         {total_mass_change_methylation:10.4f} Da")
    print(f"4. Mass of sodium for [M+Na]+ adduct:                        {mass_na:10.4f} Da")
    print("-" * 70)
    print("Final m/z Calculation:")
    print(f"m/z = (Native Glycan) + (Amidation) + (Permethylation) + (Sodium)")
    print(f"m/z = {mass_native_glycan:.4f} + ({mass_change_amidation:.4f}) + {total_mass_change_methylation:.4f} + {mass_na:.4f}")
    print(f"Final predicted m/z for [M+Na]+ ion:                        {final_mz:10.4f}")
    print("-" * 70)

if __name__ == "__main__":
    main()
