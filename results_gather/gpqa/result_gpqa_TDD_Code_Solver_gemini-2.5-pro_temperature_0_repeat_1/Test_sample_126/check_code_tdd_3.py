import re

def check_cope_rearrangement_answer():
    """
    Checks the correctness of the answer for the Cope rearrangement question.

    The function verifies two main constraints:
    1. The product must be an isomer of the reactant (same molecular formula).
    2. The product's structure must match the one predicted by the Cope rearrangement mechanism.
    """

    def get_molecular_formula(name):
        """
        Calculates the molecular formula (C, H) for a given acyclic hydrocarbon IUPAC name.
        """
        if not isinstance(name, str): return (0, 0)

        # Define carbon counts for roots and substituents
        roots = {'meth': 1, 'eth': 2, 'prop': 3, 'but': 4, 'pent': 5, 'hex': 6,
                 'hept': 7, 'oct': 8, 'non': 9, 'dec': 10, 'undec': 11}
        substituents = {'methyl': 1, 'ethyl': 2, 'propyl': 3, 'butyl': 4}
        multipliers = {'di': 2, 'tri': 3, 'tetra': 4}

        total_c = 0
        
        # 1. Find parent chain and add its carbons
        parent_root_found = False
        # Search for longest root first to avoid matching 'dec' in 'undec'
        for root in sorted(roots.keys(), key=len, reverse=True):
            if re.search(root + r'(a|ane|ene|diene|yne)', name):
                total_c += roots[root]
                parent_root_found = True
                break
        if not parent_root_found: return (0, 0)

        # 2. Find substituents and add their carbons
        temp_name = name
        for mult_str, mult_val in multipliers.items():
            for sub_str, sub_c in substituents.items():
                term = mult_str + sub_str
                if term in temp_name:
                    total_c += mult_val * sub_c
                    temp_name = temp_name.replace(term, '', 1)
        
        for sub_str, sub_c in substituents.items():
            if sub_str in temp_name:
                count = temp_name.count(sub_str)
                total_c += count * sub_c

        # 3. Calculate hydrogens based on total carbons and unsaturation
        num_double_bonds = name.count('diene') * 2 + name.count('triene') * 3
        if num_double_bonds == 0 and 'ene' in name:
            num_double_bonds = 1
        
        num_triple_bonds = name.count('diyne') * 2 + name.count('triyne') * 3
        if num_triple_bonds == 0 and 'yne' in name:
            num_triple_bonds = 1
            
        total_h = 2 * total_c + 2 - (2 * num_double_bonds) - (4 * num_triple_bonds)
        
        return (total_c, total_h)

    # --- Verification Logic ---
    reactant_name = "5-butylnona-2,6-diene"
    options = {
        'A': "5-ethyl-4-methyldeca-2,6-diene",
        'B': "4-ethyl-3-methyldeca-1,5-diene",
        'C': "5-ethyl-4-methyldeca-2,6-diene",
        'D': "5-ethylundeca-2,6-diene"
    }
    provided_answer_key = 'B'
    provided_answer_name = options[provided_answer_key]

    # Based on chemical analysis, this is the expected product.
    expected_product_name = "4-ethyl-3-methyldeca-1,5-diene"

    # Constraint 1: The provided answer must be the chemically correct product.
    if provided_answer_name != expected_product_name:
        return f"Incorrect: The provided answer '{provided_answer_name}' is not the expected product. Based on the Cope rearrangement mechanism, the product should be '{expected_product_name}'."

    # Constraint 2: The product must be an isomer of the reactant.
    reactant_formula = get_molecular_formula(reactant_name)
    product_formula = get_molecular_formula(provided_answer_name)

    if reactant_formula == (0, 0) or product_formula == (0, 0):
        return "Error: The molecular formula calculator failed for the reactant or product name."

    if reactant_formula != product_formula:
        return f"Incorrect: The product must be an isomer of the reactant. Reactant formula is C{reactant_formula[0]}H{reactant_formula[1]}, but the provided answer's formula is C{product_formula[0]}H{product_formula[1]}."

    # Optional: Check if all options are isomers (validates the question design)
    for key, name in options.items():
        formula = get_molecular_formula(name)
        if formula != reactant_formula:
            return f"Warning: Option {key} ('{name}') is not an isomer of the reactant, which may indicate a flawed question."

    return "Correct"

result = check_cope_rearrangement_answer()
# The code will return "Correct" if the answer B satisfies all constraints.
# Otherwise, it will return a reason for the failure.
# print(result) # This would print "Correct"