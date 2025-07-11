import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a sodiated, amidated, and permethylated A2G2S2 glycan.
    """
    # Step 1: Define monoisotopic atomic masses
    atomic_mass = {
        'H': 1.007825,
        'C': 12.000000,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770,
    }

    # Define formulas for monosaccharide residues
    residue_formula = {
        'Man': collections.Counter({'C': 6, 'H': 10, 'O': 5}), # Mannose
        'Gal': collections.Counter({'C': 6, 'H': 10, 'O': 5}), # Galactose
        'GlcNAc': collections.Counter({'C': 8, 'H': 13, 'N': 1, 'O': 5}), # N-acetylglucosamine
        'Neu5Ac': collections.Counter({'C': 11, 'H': 17, 'N': 1, 'O': 8}), # N-acetylneuraminic acid (Sialic Acid)
        'H2O': collections.Counter({'H': 2, 'O': 1}),
    }

    def get_mass(formula_counter):
        """Calculates mass from a formula represented by a Counter."""
        mass = 0
        for atom, count in formula_counter.items():
            mass += count * atomic_mass[atom]
        return mass

    # Step 2: Calculate the mass of the native A2G2S2 glycan
    # Composition: 3 Man + 4 GlcNAc + 2 Gal + 2 Neu5Ac
    native_glycan_formula = (3 * residue_formula['Man'] + 
                             4 * residue_formula['GlcNAc'] + 
                             2 * residue_formula['Gal'] + 
                             2 * residue_formula['Neu5Ac'] + 
                             residue_formula['H2O']) # Add water for reducing end
    
    native_glycan_mass = get_mass(native_glycan_formula)

    # Step 3: Calculate mass change from amidation
    # Reaction: R-COOH -> R-CONH2. Net change is replacing -OH with -NH2.
    # This happens on the two sialic acids.
    mass_change_amidation = get_mass(collections.Counter({'N': 1, 'H': 2})) - get_mass(collections.Counter({'O': 1, 'H': 1}))
    mass_after_amidation = native_glycan_mass + 2 * mass_change_amidation

    # Step 4: Calculate mass gain from permethylation
    # Count methylation sites on the amidated glycan:
    # - Man/Gal (non-terminal): 2 -OH
    # - GlcNAc (non-terminal): 2 -OH, 1 -NH
    # - For A2G2S2, a full count is complex. We use a known standard count and adjust.
    #   A native A2G2S2 has 38 methylation sites (OH, NH, and COOH->COOCH3).
    #   Our reaction changes COOH to CONH2.
    #   -COOH -> -COOCH3 is 1 methylation site.
    #   -CONH2 -> -CON(CH3)2 is 2 methylation sites.
    #   So, for each of the two sialic acids, we gain 1 methylation site.
    #   Total sites = 38 (native) + 2 (from amidation) = 40 sites.
    num_methylation_sites = 40
    
    # Mass gain per site is one methyl group (-CH3) replacing a proton (-H). Net is CH2.
    mass_gain_per_methylation = get_mass(collections.Counter({'C': 1, 'H': 2}))
    total_mass_gain_methylation = num_methylation_sites * mass_gain_per_methylation
    
    final_neutral_mass = mass_after_amidation + total_mass_gain_methylation

    # Step 5: Calculate the m/z of the singly-sodiated ion [M+Na]+
    sodiuated_ion_mass = final_neutral_mass + atomic_mass['Na']

    # --- Output Results ---
    print("The expected m/z for the modified glycans is calculated as follows:")
    print("This mass is identical for all three isomers: A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2.\n")
    print("Calculation Breakdown:")
    print(f"1. Mass of Native Glycan (A2G2S2): {native_glycan_mass:.4f} Da")
    print(f"2. Mass Change from Amidation (2 x -OH -> -NH2): 2 * {mass_change_amidation:.4f} Da = {2*mass_change_amidation:.4f} Da")
    print(f"3. Mass Gain from Permethylation ({num_methylation_sites} sites): {num_methylation_sites} * {mass_gain_per_methylation:.4f} Da = {total_mass_gain_methylation:.4f} Da")
    print(f"4. Mass of Sodium Ion (Na+): {atomic_mass['Na']:.4f} Da\n")
    
    print("Final Equation:")
    print(f"{native_glycan_mass:.4f} + {2*mass_change_amidation:.4f} + {total_mass_gain_methylation:.4f} + {atomic_mass['Na']:.4f} = {sodiuated_ion_mass:.4f}\n")
    
    print("Therefore, the expected m/z for the [M+Na]+ ion for each glycan is:")
    print(f"A2G(4)2S(3)2: {sodiuated_ion_mass:.4f}")
    print(f"A2G(4)S(3)S(6): {sodiuated_ion_mass:.4f}")
    print(f"A2G(4)2S(6)2: {sodiuated_ion_mass:.4f}")

calculate_glycan_mass()