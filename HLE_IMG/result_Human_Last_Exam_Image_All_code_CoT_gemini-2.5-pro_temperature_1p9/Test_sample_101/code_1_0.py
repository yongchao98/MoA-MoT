def calculate_molar_mass():
    """
    Calculates and prints the molar mass of Compound A (C13H11N3O).
    The structure of Compound A is (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile.
    """
    
    # Standard atomic weights
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }
    
    # Atom counts in Compound A (C13H11N3O)
    atom_counts = {
        'C': 13,
        'H': 11,
        'N': 3,
        'O': 1
    }
    
    # Calculate molar mass
    molar_mass = (atom_counts['C'] * atomic_weights['C'] +
                   atom_counts['H'] * atomic_weights['H'] +
                   atom_counts['N'] * atomic_weights['N'] +
                   atom_counts['O'] * atomic_weights['O'])

    # Print the detailed calculation equation
    c_mass = atom_counts['C'] * atomic_weights['C']
    h_mass = atom_counts['H'] * atomic_weights['H']
    n_mass = atom_counts['N'] * atomic_weights['N']
    o_mass = atom_counts['O'] * atomic_weights['O']
    
    print("The final product, Compound A, is (phenylamino)(3-hydroxypyridin-2-yl)acetonitrile.")
    print("Its chemical formula is C13H11N3O.")
    print("\nCalculating the molar mass of Compound A:")
    print(f"Molar Mass = ({atom_counts['C']} * {atomic_weights['C']}) + ({atom_counts['H']} * {atomic_weights['H']}) + ({atom_counts['N']} * {atomic_weights['N']}) + ({atom_counts['O']} * {atomic_weights['O']})")
    print(f"Molar Mass = {c_mass:.3f} + {h_mass:.3f} + {n_mass:.3f} + {o_mass:.3f}")
    print(f"Molar Mass = {molar_mass:.3f} g/mol")

calculate_molar_mass()