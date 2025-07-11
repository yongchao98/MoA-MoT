def solve_reaction():
    """
    This script determines the properties of Compound A from the described reaction.

    The reaction is a two-step process:
    1. Imine formation: 3-hydroxy-pyridine-2-carbaldehyde + aniline -> imine intermediate.
    2. Nucleophilic addition: The imine intermediate reacts with NaCN in a Strecker-type reaction
       to form an alpha-aminonitrile, which is Compound A.

    Compound A is 2-((cyano)(phenylamino)methyl)pyridin-3-ol.
    This script calculates its molecular formula and molar mass.
    """

    # Atomic masses (g/mol)
    atomic_masses = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Atom count in Compound A (2-((cyano)(phenylamino)methyl)pyridin-3-ol)
    # Pyridine ring part (C5H3N) + OH group = C5H4NO
    # Phenylamino group (C6H5NH) = C6H6N
    # Methine-cyano part (-CH(CN)-) = C2H
    # Let's count atoms from the final structure:
    # Carbons: 5 (pyridine ring) + 1 (methine C) + 6 (phenyl ring) + 1 (cyano C) = 13
    # Hydrogens: 3 (pyridine ring) + 1 (OH) + 1 (methine H) + 5 (phenyl ring) + 1 (NH) = 11
    # Nitrogens: 1 (pyridine ring) + 1 (amino N) + 1 (cyano N) = 3
    # Oxygens: 1 (OH) = 1
    atom_counts = {
        'C': 13,
        'H': 11,
        'N': 3,
        'O': 1
    }

    # Calculate molar mass
    molar_mass = 0
    molecular_formula = ""
    for element, count in atom_counts.items():
        molar_mass += count * atomic_masses[element]
        molecular_formula += f"{element}{count}"

    print("Compound A is 2-((cyano)(phenylamino)methyl)pyridin-3-ol.")
    print(f"The molecular formula of Compound A is: {molecular_formula}")
    print(f"The molar mass of Compound A is: {molar_mass:.3f} g/mol")

solve_reaction()