def calculate_glycan_mz():
    """
    Calculates the m/z of a permethylated, amidated, sodiated A2G2S2 N-glycan.
    """
    # Monoisotopic atomic masses
    ATOMIC_MASS = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # 1. Calculate the mass of the base, underivatized glycan (4xGlcNAc, 3xMan, 2xGal, 2xNeu5Ac)
    # The molecular formula for an A2G2S2 N-glycan (C84 H138 N6 O62) is calculated from its
    # constituent monosaccharides minus the H2O molecules lost forming 10 glycosidic bonds.
    glycan_C = 84
    glycan_H = 138
    glycan_N = 6
    glycan_O = 62
    base_glycan_mass = (glycan_C * ATOMIC_MASS['C'] +
                        glycan_H * ATOMIC_MASS['H'] +
                        glycan_N * ATOMIC_MASS['N'] +
                        glycan_O * ATOMIC_MASS['O'])

    # 2. Calculate the mass change from amidation of two sialic acids
    # The reaction -COOH -> -CONH2 corresponds to a net change of -O +NH per site.
    # Since there are two sialic acids, we multiply the change by 2.
    amidation_mass_change = 2 * (ATOMIC_MASS['N'] + ATOMIC_MASS['H'] - ATOMIC_MASS['O'])

    # 3. Calculate the mass increase from permethylation
    # We count the number of sites (-OH and -NH2) available for methylation.
    # - Monomer sites: GlcNAc (4 OH), Hexose (5 OH), Neu5Ac-amide (5 OH + 1 NH2 -> 7 sites)
    # - Total sites on monomers = 4*4 (GlcNAc) + 5*5 (Hexoses) + 2*7 (Neu5Ac-amides) = 16 + 25 + 14 = 55 sites
    # - Number of linkages = 11 residues - 1 = 10. Each linkage removes 2 sites.
    # - Total methylations = 55 - (2 * 10) = 35 sites.
    # The mass increase per methylation (replacing -H with -CH3) is the mass of -CH2.
    num_methylations = 35
    permethylation_mass_increase = num_methylations * (ATOMIC_MASS['C'] + 2 * ATOMIC_MASS['H'])
    
    # 4. Sum all components and add the mass of the sodium adduct
    sodium_mass = ATOMIC_MASS['Na']
    final_mz = base_glycan_mass + amidation_mass_change + permethylation_mass_increase + sodium_mass

    # Output the final equation as requested
    print(f"The m/z for the glycans is calculated as follows:")
    print(f"{base_glycan_mass:.4f} (Base Glycan) + {amidation_mass_change:.4f} (Amidation) + {permethylation_mass_increase:.4f} (Permethylation) + {sodium_mass:.4f} (Sodium Adduct) = {final_mz:.4f} m/z")

calculate_glycan_mz()