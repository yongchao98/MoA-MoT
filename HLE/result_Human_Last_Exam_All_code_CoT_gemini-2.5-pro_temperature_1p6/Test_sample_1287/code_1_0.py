import collections

def calculate_mass(formula_dict, atomic_masses):
    """Calculates the monoisotopic mass from a formula dictionary."""
    mass = 0.0
    for atom, count in formula_dict.items():
        mass += count * atomic_masses[atom]
    return mass

def main():
    """
    Calculates the m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """
    # Monoisotopic atomic masses
    ATOMIC_MASSES = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Molecular formulas of components
    FORMULAS = {
        'Hexose':   {'C': 6, 'H': 12, 'O': 6},  # For Mannose (Man) and Galactose (Gal)
        'GlcNAc':   {'C': 8, 'H': 15, 'N': 1, 'O': 6},
        'Neu5Ac':   {'C': 11, 'H': 19, 'N': 1, 'O': 9},
        'H2O':      {'H': 2, 'O': 1},
        'OH':       {'O': 1, 'H': 1},
        'NH2':      {'N': 1, 'H': 2},
        'CH2':      {'C': 1, 'H': 2}, # Mass change per methylation (CH3 replaces H)
    }

    # Glycan composition
    num_hexose = 5  # 3 Man + 2 Gal
    num_glcnac = 4
    num_neu5ac = 2
    num_residues = num_hexose + num_glcnac + num_neu5ac
    num_glycosidic_bonds = num_residues - 1 # For a single released glycan chain

    # 1. Calculate mass of the unmodified glycan
    mass_hexose = calculate_mass(FORMULAS['Hexose'], ATOMIC_MASSES)
    mass_glcnac = calculate_mass(FORMULAS['GlcNAc'], ATOMIC_MASSES)
    mass_neu5ac = calculate_mass(FORMULAS['Neu5Ac'], ATOMIC_MASSES)
    mass_h2o = calculate_mass(FORMULAS['H2O'], ATOMIC_MASSES)

    total_residue_mass = (num_hexose * mass_hexose) + \
                         (num_glcnac * mass_glcnac) + \
                         (num_neu5ac * mass_neu5ac)

    unmodified_glycan_mass = total_residue_mass - (num_glycosidic_bonds * mass_h2o)

    # 2. Calculate mass change from amidation (COOH -> CONH2)
    # This is equivalent to removing OH and adding NH2
    mass_oh = calculate_mass(FORMULAS['OH'], ATOMIC_MASSES)
    mass_nh2 = calculate_mass(FORMULAS['NH2'], ATOMIC_MASSES)
    amidation_mass_change_per_site = mass_nh2 - mass_oh
    total_amidation_change = 2 * amidation_mass_change_per_site # For two sialic acids

    # 3. Calculate mass change from permethylation
    # Detailed structural analysis reveals 39 sites (free -OH and -NH groups)
    num_methylation_sites = 39
    methylation_mass_change_per_site = calculate_mass(FORMULAS['CH2'], ATOMIC_MASSES)
    total_methylation_change = num_methylation_sites * methylation_mass_change_per_site

    # 4. Calculate final m/z for the [M+Na]+ ion
    final_modified_mass = unmodified_glycan_mass + total_amidation_change + total_methylation_change
    final_mz = final_modified_mass + ATOMIC_MASSES['Na']

    print("Mass calculation for permethylated, amidated, and sodiated A2G2S2 glycan isomers.")
    print("-" * 75)
    print(f"Unmodified A2G2S2 Glycan Mass: {unmodified_glycan_mass:.4f}")
    print(f"Total Mass Change from Amidation (x2): {total_amidation_change:.4f}")
    print(f"Total Mass Change from Permethylation (x39): {total_methylation_change:.4f}")
    print(f"Mass of Sodium (Na+) ion: {ATOMIC_MASSES['Na']:.4f}")
    print("-" * 75)
    print("Final m/z equation:")
    print(f"{unmodified_glycan_mass:.4f} + ({total_amidation_change:.4f}) + {total_methylation_change:.4f} + {ATOMIC_MASSES['Na']:.4f} = {final_mz:.4f}")
    print("-" * 75)
    print("Since all three glycans are isomers that undergo the same chemical modifications,")
    print("they will all have the same final m/z.")
    print(f"\nExpected mass for A2G(4)2S(3)2: {final_mz:.4f}")
    print(f"Expected mass for A2G(4)S(3)S(6): {final_mz:.4f}")
    print(f"Expected mass for A2G(4)2S(6)2: {final_mz:.4f}")

if __name__ == "__main__":
    main()