import collections

def calculate_mass(formula, atomic_masses):
    """Calculates the monoisotopic mass from an elemental formula dictionary."""
    mass = 0.0
    for element, count in formula.items():
        mass += count * atomic_masses[element]
    return mass

def main():
    # Monoisotopic masses of elements and common groups
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915,
        'Na': 22.989770
    }
    
    # Elemental composition of the starting A2G2S2 glycan: C84 H138 N6 O62
    base_glycan_formula = collections.defaultdict(int, {'C': 84, 'H': 138, 'N': 6, 'O': 62})
    base_glycan_mass = calculate_mass(base_glycan_formula, atomic_masses)

    # Initial number of methylatable sites (OH and NH groups) in the base glycan
    initial_methylation_sites = 31

    # --- Glycan 1: A2G(4)2S(3)2 (both sialic acids are alpha-2,3 linked) ---
    # Both sialic acids undergo lactamization (-H2O each)
    glycan1_formula = base_glycan_formula.copy()
    glycan1_formula['H'] -= 4  # 2 * H2
    glycan1_formula['O'] -= 2  # 2 * O
    
    # Each lactamization consumes one methylation site (the C8-OH)
    glycan1_methylation_sites = initial_methylation_sites - 2
    
    # Permethylation adds CH2 per site (replaces H with CH3)
    glycan1_formula['C'] += glycan1_methylation_sites
    glycan1_formula['H'] += glycan1_methylation_sites * 2
    
    # Calculate final mass
    neutral_mass_1 = calculate_mass(glycan1_formula, atomic_masses)
    sodiated_mass_1 = neutral_mass_1 + atomic_masses['Na']

    print("--- Glycan 1: A2G(4)2S(3)2 (2x α2,3) ---")
    print(f"This glycan undergoes two lactamization reactions and has {glycan1_methylation_sites} permethylation sites.")
    print(f"Final sodiated mass [M+Na]+ = Neutral Mass ({neutral_mass_1:.4f} Da) + Na Mass ({atomic_masses['Na']:.4f} Da) = {sodiated_mass_1:.4f} m/z\n")

    # --- Glycan 2: A2G(4)S(3)S(6) (one alpha-2,3 and one alpha-2,6) ---
    # One lactamization (-H2O) and one amidation (-O +NH)
    glycan2_formula = base_glycan_formula.copy()
    glycan2_formula['H'] -= 2  # from lactamization
    glycan2_formula['O'] -= 1  # from lactamization
    glycan2_formula['O'] -= 1  # from amidation
    glycan2_formula['N'] += 1  # from amidation
    glycan2_formula['H'] += 1  # from amidation
    
    # Lactamization consumes one site, amidation creates two sites (on -CONH2)
    glycan2_methylation_sites = initial_methylation_sites - 1 + 2
    
    # Permethylation adds CH2 per site
    glycan2_formula['C'] += glycan2_methylation_sites
    glycan2_formula['H'] += glycan2_methylation_sites * 2

    # Calculate final mass
    neutral_mass_2 = calculate_mass(glycan2_formula, atomic_masses)
    sodiated_mass_2 = neutral_mass_2 + atomic_masses['Na']

    print("--- Glycan 2: A2G(4)S(3)S(6) (1x α2,3, 1x α2,6) ---")
    print(f"This glycan undergoes one lactamization and one amidation reaction, and has {glycan2_methylation_sites} permethylation sites.")
    print(f"Final sodiated mass [M+Na]+ = Neutral Mass ({neutral_mass_2:.4f} Da) + Na Mass ({atomic_masses['Na']:.4f} Da) = {sodiated_mass_2:.4f} m/z\n")

    # --- Glycan 3: A2G(4)2S(6)2 (both sialic acids are alpha-2,6 linked) ---
    # Both sialic acids undergo amidation (-O +NH each)
    glycan3_formula = base_glycan_formula.copy()
    glycan3_formula['O'] -= 2  # 2 * O
    glycan3_formula['N'] += 2  # 2 * N
    glycan3_formula['H'] += 2  # 2 * H
    
    # Each amidation creates two methylation sites
    glycan3_methylation_sites = initial_methylation_sites + 2 + 2
    
    # Permethylation adds CH2 per site
    glycan3_formula['C'] += glycan3_methylation_sites
    glycan3_formula['H'] += glycan3_methylation_sites * 2
    
    # Calculate final mass
    neutral_mass_3 = calculate_mass(glycan3_formula, atomic_masses)
    sodiated_mass_3 = neutral_mass_3 + atomic_masses['Na']
    
    print("--- Glycan 3: A2G(4)2S(6)2 (2x α2,6) ---")
    print(f"This glycan undergoes two amidation reactions and has {glycan3_methylation_sites} permethylation sites.")
    print(f"Final sodiated mass [M+Na]+ = Neutral Mass ({neutral_mass_3:.4f} Da) + Na Mass ({atomic_masses['Na']:.4f} Da) = {sodiated_mass_3:.4f} m/z")


if __name__ == '__main__':
    main()
<<<2616.2055, 2675.2790, 2734.3525>>>