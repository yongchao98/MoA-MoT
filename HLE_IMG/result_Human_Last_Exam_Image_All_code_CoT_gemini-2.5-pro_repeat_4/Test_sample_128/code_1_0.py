def solve_chemistry_problem():
    """
    Identifies compound A from the reaction and calculates its molar mass.
    """
    # Product Information
    compound_name = "1-cyano-2-(pyridin-2-yl)isoindole"
    # The molecular formula is derived from the structure of the product:
    # It consists of an isoindole core (C8H6N), a cyano group (CN) replacing an H,
    # and a pyridin-2-yl group (C5H4N) replacing the H on the nitrogen.
    # Total: C(8+1+5) H(6-1-1+4) N(1+1+1) = C14H8N3 -> Wait, let's recount H.
    # Isoindole (C8H7N). H on C1 replaced by CN. H on N replaced by Py.
    # C8H5N + CN + C5H4N = C(8+1+5) H(5+4) N(1+1+1) = C14H9N3.
    molecular_formula = "C14H9N3"
    
    # Atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
    }

    # Counts of each atom in the molecule
    atom_counts = {
        'C': 14,
        'H': 9,
        'N': 3,
    }

    # Calculate molar mass
    molar_mass = (atom_counts['C'] * atomic_weights['C'] +
                   atom_counts['H'] * atomic_weights['H'] +
                   atom_counts['N'] * atomic_weights['N'])

    print(f"Compound A is: {compound_name}")
    print(f"Its molecular formula is: {molecular_formula}")
    print("\nCalculating the molar mass (g/mol):")
    
    c_mass = atom_counts['C'] * atomic_weights['C']
    h_mass = atom_counts['H'] * atomic_weights['H']
    n_mass = atom_counts['N'] * atomic_weights['N']
    
    # As requested, output each number in the final equation
    print(f"({atom_counts['C']} * {atomic_weights['C']}) + ({atom_counts['H']} * {atomic_weights['H']}) + ({atom_counts['N']} * {atomic_weights['N']})")
    print(f"= {c_mass:.3f} + {h_mass:.3f} + {n_mass:.3f}")
    print(f"= {molar_mass:.3f} g/mol")

solve_chemistry_problem()