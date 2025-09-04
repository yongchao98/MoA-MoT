import re

def check_answer():
    """
    Checks the correctness of the given answer for the chemistry question.
    The logic involves two main checks:
    1. Isomer Check: The product must have the same molecular formula as the reactant.
    2. Reaction Mechanism Check: The product must be the one predicted by the Cope rearrangement rules.
    """
    
    # --- Problem Definition ---
    reactant_name = "5-butylnona-2,6-diene"
    llm_answer_choice = "B"
    options = {
        "A": "5-ethyl-4-methyldeca-2,6-diene",
        "B": "4-ethyl-3-methyldeca-1,5-diene",
        "C": "5-ethyl-4-methyldeca-2,6-diene",
        "D": "5-ethylundeca-2,6-diene"
    }
    llm_answer_name = options[llm_answer_choice]

    # --- Helper Functions ---

    def get_molecular_formula(name):
        """
        Calculates the molecular formula (C, H) for simple acyclic hydrocarbons from their IUPAC name.
        """
        if not isinstance(name, str):
            return None
            
        parents = {'undec': 11, 'dec': 10, 'non': 9, 'oct': 8, 'hept': 7, 'hex': 6, 'pent': 5, 'but': 4, 'prop': 3, 'eth': 2, 'meth': 1}
        substituents = {'butyl': 4, 'propyl': 3, 'ethyl': 2, 'methyl': 1}
        multipliers = {'di': 2, 'tri': 3, 'tetra': 4}

        total_c = 0
        
        # 1. Find and add carbons from the parent chain
        parent_found = False
        # Sort by length to find 'undec' before 'dec', etc.
        for p_name, p_c in sorted(parents.items(), key=lambda item: len(item[0]), reverse=True):
            # Parent chain is usually at the end of a segment
            if re.search(f'{p_name}(a|ane|ene|diene|yne|triene)?$', name):
                total_c += p_c
                parent_found = True
                break
        if not parent_found:
            return None # Could not parse parent chain

        # 2. Find and add carbons from substituents
        temp_name = name
        for m_name, m_val in multipliers.items():
            for s_name, s_c in substituents.items():
                term = m_name + s_name
                if term in temp_name:
                    total_c += m_val * s_c
                    temp_name = temp_name.replace(term, '', 1)
        
        for s_name, s_c in substituents.items():
            count = temp_name.count(s_name)
            if count > 0:
                total_c += count * s_c

        # 3. Calculate hydrogens based on total carbons and degrees of unsaturation
        total_h = 2 * total_c + 2  # Start with alkane formula
        
        if 'triene' in name:
            total_h -= 6
        elif 'diene' in name:
            total_h -= 4
        elif 'yne' in name: # Triple bond = 2 pi bonds
            total_h -= 4
        elif 'ene' in name: # Double bond = 1 pi bond
            total_h -= 2
            
        return (total_c, total_h)

    def predict_cope_product(reactant):
        """
        Predicts the product of a Cope rearrangement for a specific, known reactant.
        This function hardcodes the result of the chemical analysis.
        """
        if reactant == "5-butylnona-2,6-diene":
            # Based on the mechanism: C4-C5 bond breaks, C2-C7 bond forms, pi bonds shift.
            # This results in a new 10-carbon chain (deca-1,5-diene) with
            # an ethyl group at C4 and a methyl group at C3.
            return "4-ethyl-3-methyldeca-1,5-diene"
        return None

    # --- Verification Logic ---

    # 1. Isomer Check
    reactant_formula = get_molecular_formula(reactant_name)
    answer_formula = get_molecular_formula(llm_answer_name)

    if reactant_formula is None:
        return "Checker Error: Could not determine molecular formula for the reactant."
    if answer_formula is None:
        return f"Checker Error: Could not determine molecular formula for the answer '{llm_answer_name}'."

    if reactant_formula != answer_formula:
        return (f"Incorrect. The product of a rearrangement must be an isomer of the reactant. "
                f"Reactant formula is C{reactant_formula[0]}H{reactant_formula[1]}, but the "
                f"answer's formula is C{answer_formula[0]}H{answer_formula[1]}.")

    # 2. Reaction Mechanism Check
    correct_product = predict_cope_product(reactant_name)
    
    if correct_product is None:
        return "Checker Error: The reaction logic for this reactant is not implemented."

    if llm_answer_name == correct_product:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer '{llm_answer_name}' is not the correct product. "
                f"The Cope rearrangement of {reactant_name} yields '{correct_product}'.")

# Run the check and print the result
result = check_answer()
print(result)