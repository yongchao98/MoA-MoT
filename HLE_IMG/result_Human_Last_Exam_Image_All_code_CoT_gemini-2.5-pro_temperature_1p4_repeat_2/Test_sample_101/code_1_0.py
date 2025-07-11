def calculate_molecular_formula(reactant1_formula, reactant2_formula, loss_formula, addition_formula):
    """
    Calculates the molecular formula of the final product from a two-step reaction.
    Step 1: Condensation of two reactants with a loss of a small molecule.
    Step 2: Addition of a group to the intermediate.
    """
    # Helper function to parse a formula string into a dictionary
    def parse_formula(formula):
        import re
        atom_counts = {}
        # Find all element-count pairs (e.g., C6, H5, N, O2)
        # N is treated as N1
        for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
            atom_counts[element] = atom_counts.get(element, 0) + (int(count) if count else 1)
        return atom_counts

    # Helper function to combine or subtract atom dictionaries
    def combine_atoms(d1, d2, operation='+'):
        result = d1.copy()
        for element, count in d2.items():
            if operation == '+':
                result[element] = result.get(element, 0) + count
            elif operation == '-':
                result[element] = result.get(element, 0) - count
        return {k: v for k, v in result.items() if v > 0} # Remove elements with zero or negative count

    # Helper function to format a dictionary back into a formula string
    def format_formula(atom_dict):
        # Order elements C, H, then alphabetically
        elements = sorted(atom_dict.keys(), key=lambda x: ('0' if x == 'C' else ('1' if x == 'H' else x)))
        formula_str = ""
        for element in elements:
            count = atom_dict[element]
            formula_str += element + (str(count) if count > 1 else "")
        return formula_str

    # Parse initial formulas
    reactant1 = parse_formula(reactant1_formula)
    reactant2 = parse_formula(reactant2_formula)
    molecule_lost = parse_formula(loss_formula)
    molecule_added = parse_formula(addition_formula)

    print(f"Reactant 1 (3-hydroxy-pyridine-2-carbaldehyde): {format_formula(reactant1)}")
    print(f"Reactant 2 (Aniline): {format_formula(reactant2)}")
    print("-" * 30)
    
    # Step 1: Imine formation
    # Combine reactants and subtract water
    print("Step 1: Imine Formation (Condensation)")
    combined_reactants = combine_atoms(reactant1, reactant2, operation='+')
    print(f"Combined Reactants: {format_formula(combined_reactants)}")
    print(f"Molecule Lost: {format_formula(molecule_lost)}")
    intermediate = combine_atoms(combined_reactants, molecule_lost, operation='-')
    print(f"Intermediate Formula: {format_formula(intermediate)}")
    print("-" * 30)

    # Step 2: Cyanide addition
    # Add HCN to the intermediate
    print("Step 2: Cyanide Addition (Strecker-type)")
    print(f"Group Added (from NaCN + proton source): {format_formula(molecule_added)}")
    compound_a = combine_atoms(intermediate, molecule_added, operation='+')
    
    print("-" * 30)
    print("Final Equation for Compound A:")
    # Print the final counts explicitly as requested
    final_formula_dict = compound_a
    final_formula_string = ""
    elements = sorted(final_formula_dict.keys(), key=lambda x: ('0' if x == 'C' else ('1' if x == 'H' else x)))
    for element in elements:
        final_formula_string += f"{element}{final_formula_dict[element]} "
        
    print(f"The molecular formula of Compound A is: {final_formula_string.strip()}")


# --- Main execution ---
# Define molecular formulas for reactants and reaction steps
r1_formula = "C6H5NO2"   # 3-hydroxy-pyridine-2-carbaldehyde
r2_formula = "C6H7N"     # Aniline
h2o_formula = "H2O"      # Water lost during condensation
hcn_formula = "HCN"      # Effective group added from NaCN + proton source

calculate_molecular_formula(r1_formula, r2_formula, h2o_formula, hcn_formula)