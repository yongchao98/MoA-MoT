import collections

def get_formula(atom_counts):
    """Converts a dictionary of atom counts into a chemical formula string."""
    # Standard order for organic compounds: C, H, then alphabetical for the rest
    formula = ""
    if 'C' in atom_counts:
        count = atom_counts['C']
        formula += 'C'
        if count > 1:
            formula += str(count)
    if 'H' in atom_counts:
        count = atom_counts['H']
        formula += 'H'
        if count > 1:
            formula += str(count)
    
    # Add other elements alphabetically
    other_elements = sorted([elem for elem in atom_counts if elem not in ['C', 'H']])
    for elem in other_elements:
        count = atom_counts[elem]
        formula += elem
        if count > 1:
            formula += str(count)
    return formula

def main():
    """
    Calculates and displays the product of the reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    # Step 1: Define the reactants using dictionaries for atom counts.
    butadiene = {'C': 4, 'H': 6}
    dienophile = {'C': 2, 'Cl': 2, 'F': 2}

    # Step 2: The reaction is a 1:1 addition (Diels-Alder).
    # We use a Counter to easily add the atom counts.
    product_atoms = collections.Counter(butadiene) + collections.Counter(dienophile)

    # Step 3: Format the chemical formulas from the atom counts.
    butadiene_formula = get_formula(butadiene)
    dienophile_formula = get_formula(dienophile)
    product_formula = get_formula(product_atoms)
    
    # Step 4: Print the explanation and the final equation.
    print("The reaction between butadiene and 1,1-dichloro-2,2-difluoroethene is a Diels-Alder [4+2] cycloaddition.")
    print("This reaction forms a new six-membered ring.")
    print("\nThe final reaction equation is:")
    
    # The equation contains all the required numbers within the molecular formulas.
    # Reactant 1 + Reactant 2 -> Product
    equation = f"{butadiene_formula} + {dienophile_formula} -> {product_formula}"
    print(equation)
    
    print("\nProduct Name: 4,4-dichloro-5,5-difluorocyclohex-1-ene")
    
if __name__ == "__main__":
    main()