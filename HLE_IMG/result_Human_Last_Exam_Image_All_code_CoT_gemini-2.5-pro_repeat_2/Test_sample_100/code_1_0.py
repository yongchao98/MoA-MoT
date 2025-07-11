def add_formulas(formula1, formula2):
    """Combines two molecular formulas represented as dictionaries."""
    result = formula1.copy()
    for atom, count in formula2.items():
        result[atom] = result.get(atom, 0) + count
    return result

def format_formula(formula):
    """Formats a dictionary of atoms into a standard molecular formula string."""
    # Standard order: C, H, then alphabetical for the rest
    order = ['C', 'H']
    other_atoms = sorted([atom for atom in formula.keys() if atom not in order])
    full_order = order + other_atoms
    
    formula_str = ""
    for atom in full_order:
        if atom in formula:
            count = formula[atom]
            formula_str += atom
            if count > 1:
                formula_str += str(count)
    return formula_str

# Based on the analysis, Intermediate 1 is 5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
# Its molecular formula is calculated as follows:
# Core skeleton (C7), Hydrogens (H6), Fluorines (F3), Nitrogen (N1), Oxygen (O1).
intermediate_1_formula = {'C': 7, 'H': 6, 'F': 3, 'N': 1, 'O': 1}

# Step 2: Hydrogenation adds two hydrogen atoms (H2).
hydrogen_molecule = {'H': 2}
intermediate_2_formula = add_formulas(intermediate_1_formula, hydrogen_molecule)

# Step 3: Hydrolysis adds one water molecule (H2O).
water_molecule = {'H': 2, 'O': 1}
product_formula = add_formulas(intermediate_2_formula, water_molecule)

print("Step-by-step determination of the molecular formula:")
print(f"1. Formula after deprotection (Intermediate 1): {format_formula(intermediate_1_formula)}")
print(f"2. Formula after hydrogenation (Intermediate 2): {format_formula(intermediate_2_formula)}")
print(f"3. Formula after hydrolysis (Final Product): {format_formula(product_formula)}")
print("\n-------------------------------------------------")
print("Final Product Molecular Formula Breakdown:")

# Output each number in the final formula
C = product_formula.get('C', 0)
H = product_formula.get('H', 0)
F = product_formula.get('F', 0)
N = product_formula.get('N', 0)
O = product_formula.get('O', 0)

print(f"The final molecular formula is C{C}H{H}F{F}N{N}O{O}")
print(f"Number of Carbon (C) atoms: {C}")
print(f"Number of Hydrogen (H) atoms: {H}")
print(f"Number of Fluorine (F) atoms: {F}")
print(f"Number of Nitrogen (N) atoms: {N}")
print(f"Number of Oxygen (O) atoms: {O}")