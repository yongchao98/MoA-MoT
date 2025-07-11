def get_molecular_formula(atom_dict):
    """Calculates the molecular formula string from a dictionary of atom counts."""
    formula = ""
    for atom in ['C', 'H', 'N', 'O']:
        if atom in atom_dict and atom_dict[atom] > 0:
            formula += atom
            if atom_dict[atom] > 1:
                formula += str(atom_dict[atom])
    return formula

def solve_reaction():
    """
    Determines and presents the structures of products A, B, and C.
    """
    # Define atom counts for proposed structures based on analysis
    # Product C: N-acetyl-2-(1-pyrrolin-2-yl)pyrrolidine-2-carboxylic acid
    structure_C = {'C': 11, 'H': 16, 'N': 2, 'O': 3}

    # Product A: 7-(2-hydroxypyrrolidin-2-yl)-1-methyl-2-(methoxycarbonyl)-2,3-dihydro-1H-pyrrolo[1,2-a]pyrrole
    structure_A = {'C': 14, 'H': 20, 'N': 2, 'O': 3}

    # Product B: (Z,3S,11aS)-3-(4,5-dihydro-1H-pyrrol-2-yl)-2,3,11,11a-tetrahydro-1H,5H-pyrrolo[2,1-c][1,4]oxazepine-1,5-dione
    # This structure is derived from SM + MP - MeOH.
    structure_B = {'C': 12, 'H': 14, 'N': 2, 'O': 3}
    
    # Get formula strings
    formula_A = get_molecular_formula(structure_A)
    formula_B = get_molecular_formula(structure_B)
    formula_C = get_molecular_formula(structure_C)

    print("Based on the analysis of the reaction and spectroscopic data, the proposed structures are:")
    print("-" * 40)
    
    print("Product A:")
    print("  Proposed Structure Name: 7-(2-hydroxypyrrolidin-2-yl)-1-methyl-2-(methoxycarbonyl)-2,3-dihydro-1H-pyrrolo[1,2-a]pyrrole")
    print(f"  Molecular Formula: {formula_A}")
    print("  Formation Pathway: [N-acetyl-SM (C)] + [Methyl Propiolate] -> [A] + CO2")
    print(f"  Equation: C11H16N2O3 + C4H4O2 -> {formula_A} + CO2")
    print("-" * 40)

    print("Product B:")
    print("  Proposed Structure Name: (Z,3S,11aS)-3-(4,5-dihydro-1H-pyrrol-2-yl)-2,3,11,11a-tetrahydro-1H,5H-pyrrolo[2,1-c][1,4]oxazepine-1,5-dione")
    print(f"  Molecular Formula: {formula_B}")
    print("  Formation Pathway: [SM] + [Methyl Propiolate] -> [B] + MeOH")
    print(f"  Equation: C9H14N2O2 + C4H4O2 -> {formula_B} + CH4O")
    print("-" * 40)

    print("Product C:")
    print("  Proposed Structure Name: N-acetyl-2-(4,5-dihydro-1H-pyrrol-2-yl)pyrrolidine-2-carboxylic acid")
    print(f"  Molecular Formula: {formula_C}")
    print("  Formation Pathway: N-acetylation of Starting Material (SM)")
    print(f"  Equation: C9H14N2O2 + (Ac) - H -> {formula_C}")
    print("-" * 40)

solve_reaction()