def get_reaction_summary():
    """
    Summarizes the Diels-Alder reaction between butadiene and
    1,1-dichloro-2,2-difluoroethene.
    """

    # Define reactants by their molecular formulas
    # Butadiene: C4H6
    # 1,1-dichloro-2,2-difluoroethene: C2Cl2F2
    reactants = {
        "butadiene": {"C": 4, "H": 6},
        "dienophile": {"C": 2, "Cl": 2, "F": 2}
    }

    # Diels-Alder is an addition reaction, so sum the atoms for the product formula
    product_formula_dict = {}
    for reactant in reactants.values():
        for atom, count in reactant.items():
            product_formula_dict[atom] = product_formula_dict.get(atom, 0) + count

    # Function to format a dictionary of atoms into a standard molecular formula string
    def format_formula(atom_dict):
        # Standard order is C, H, then alphabetical for other elements
        formula_parts = []
        if 'C' in atom_dict and atom_dict['C'] > 0:
            formula_parts.append(f"C{atom_dict['C'] if atom_dict['C'] > 1 else ''}")
        if 'H' in atom_dict and atom_dict['H'] > 0:
            formula_parts.append(f"H{atom_dict['H'] if atom_dict['H'] > 1 else ''}")
        
        other_atoms = sorted([atom for atom in atom_dict if atom not in ['C', 'H']])
        for atom in other_atoms:
            formula_parts.append(f"{atom}{atom_dict[atom] if atom_dict[atom] > 1 else ''}")
            
        return "".join(formula_parts)

    butadiene_str = format_formula(reactants["butadiene"])
    dienophile_str = format_formula(reactants["dienophile"])
    product_str = format_formula(product_formula_dict)
    
    # The name of the resulting product
    product_name = "4,4-dichloro-5,5-difluorocyclohex-1-ene"

    print("Reaction Type: Diels-Alder [4+2] Cycloaddition")
    print("-" * 50)
    print("Reactants:")
    print(f"  - Diene: Butadiene ({butadiene_str})")
    print(f"  - Dienophile: 1,1-dichloro-2,2-difluoroethene ({dienophile_str})")
    print("-" * 50)
    print("Reaction Equation (Molecular Formulas):")
    # This printout includes every number from the final equation as requested
    print(f"  {butadiene_str} + {dienophile_str} -> {product_str}")
    print("-" * 50)
    print("Product Name:")
    print(f"  {product_name}")


# Execute the function to print the summary
get_reaction_summary()