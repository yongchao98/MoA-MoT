import re

def get_molecular_formula(name: str) -> tuple[int, int] or None:
    """
    Calculates the molecular formula (C, H counts) from a simplified IUPAC name.
    """
    if not isinstance(name, str):
        return None

    # Define known chemical groups and their carbon counts
    parent_chains = {'undec': 11, 'dec': 10, 'non': 9, 'oct': 8, 'hept': 7, 'hex': 6, 'pent': 5, 'but': 4, 'prop': 3, 'eth': 2, 'meth': 1}
    substituents = {'butyl': 4, 'propyl': 3, 'ethyl': 2, 'methyl': 1}
    multipliers = {'di': 2, 'tri': 3, 'tetra': 4}

    c_count = 0

    # 1. Find parent chain carbon count (match longest first, e.g., 'undec' before 'dec')
    temp_name = name
    for p_name in sorted(parent_chains.keys(), key=len, reverse=True):
        if p_name in temp_name:
            c_count += parent_chains[p_name]
            # Remove parent name to avoid it being counted as a substituent (e.g., 'meth' in 'methyl')
            temp_name = temp_name.replace(p_name, '')
            break
    
    if c_count == 0:
        # Could not determine parent chain
        return None

    # 2. Find substituent carbon count
    sub_pattern = r'(di|tri|tetra)?(butyl|propyl|ethyl|methyl)'
    matches = re.findall(sub_pattern, temp_name)
    for prefix, sub_name in matches:
        multiplier = multipliers.get(prefix, 1)
        c_count += multiplier * substituents[sub_name]

    # 3. Calculate hydrogen count based on total carbons and unsaturation
    # Start with the alkane formula for the total carbon count
    h_count = 2 * c_count + 2

    # Subtract hydrogens for double or triple bonds
    num_double_bonds = 0
    if 'triene' in name:
        num_double_bonds = 3
    elif 'diene' in name:
        num_double_bonds = 2
    elif 'ene' in name:
        num_double_bonds = 1

    num_triple_bonds = 0
    if 'triyne' in name:
        num_triple_bonds = 3
    elif 'diyne' in name:
        num_triple_bonds = 2
    elif 'yne' in name:
        num_triple_bonds = 1
        
    h_count -= (2 * num_double_bonds) + (4 * num_triple_bonds)

    return c_count, h_count

def check_correctness():
    """
    Checks if the provided LLM answer is correct for the given chemistry question.
    """
    # --- Problem Definition ---
    reactant_name = "5-butylnona-2,6-diene"
    options = {
        'A': "5-ethyl-4-methyldeca-2,6-diene",
        'B': "4-ethyl-3-methyldeca-1,5-diene",
        'C': "5-ethyl-4-methyldeca-2,6-diene",
        'D': "5-ethylundeca-2,6-diene"
    }
    llm_answer_key = 'B'
    
    # --- Ground Truth ---
    # The heating of a 1,5-diene causes a Cope rearrangement.
    # The expected product is 4-ethyl-3-methyldeca-1,5-diene.
    expected_product_name = "4-ethyl-3-methyldeca-1,5-diene"

    # --- Verification ---
    
    # Check 1: Is the selected option key valid?
    if llm_answer_key not in options:
        return f"Incorrect: The answer key '{llm_answer_key}' is not one of the available options (A, B, C, D)."
    
    selected_product_name = options[llm_answer_key]

    # Check 2: Is the product an isomer of the reactant? (Necessary condition)
    reactant_formula = get_molecular_formula(reactant_name)
    product_formula = get_molecular_formula(selected_product_name)

    if reactant_formula is None:
        return f"Error: Could not parse the reactant name '{reactant_name}' to determine its formula."
    if product_formula is None:
        return f"Error: Could not parse the product name '{selected_product_name}' to determine its formula."

    if reactant_formula != product_formula:
        return (f"Incorrect: The product must be an isomer of the reactant. "
                f"Reactant ({reactant_name}) formula is C{reactant_formula[0]}H{reactant_formula[1]}, "
                f"but product ({selected_product_name}) formula is C{product_formula[0]}H{product_formula[1]}.")

    # Check 3: Is the product the correct one from the reaction mechanism?
    if selected_product_name != expected_product_name:
        correct_key = next((key for key, value in options.items() if value == expected_product_name), None)
        return (f"Incorrect: The molecule in answer {llm_answer_key} ({selected_product_name}) is not the correct product. "
                f"The Cope rearrangement of {reactant_name} yields {expected_product_name} (Option {correct_key}).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)
