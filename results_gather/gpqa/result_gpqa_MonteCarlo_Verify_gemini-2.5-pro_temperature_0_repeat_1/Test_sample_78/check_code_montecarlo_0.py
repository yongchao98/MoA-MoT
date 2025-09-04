def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer by verifying the chemical constraints.
    1. It defines the properties of the product based on the NMR data.
    2. It defines the properties of each possible starting material (the options).
    3. It checks which option is a valid precursor for the product.
    """

    # --- Step 1: Analyze Product from NMR Data ---
    # The NMR data strongly indicates a specific structure for the product.
    # 1H NMR: δ 7.08 (2H, d) & 7.71 (2H, d) -> para-substituted aromatic ring.
    # 1H NMR: Two methyl singlets (δ 2.28, 2.31) -> One is likely on the para-ring, one is part of a methyl ketone.
    # 13C NMR: δ 197.7 -> Ketone C=O.
    # Conclusion: The product contains a p-tolyl (4-methylphenyl) group.
    product_required_group = 'p-tolyl'
    
    # The starting material's formula is given.
    required_formula = 'C11H12O'

    # --- Step 2: Define the properties of the candidate compounds ---
    options = {
        'A': {
            'name': '2-methyl-3-styryloxirane',
            'formula': 'C11H12O',  # C6H5-CH=CH-CH(O)CH-CH3
            'aromatic_group': 'phenyl'
        },
        'B': {
            'name': '2-(1-phenylprop-1-en-2-yl)oxirane',
            'formula': 'C11H12O',  # C6H5-CH=C(CH3)-CH(O)CH2
            'aromatic_group': 'phenyl'
        },
        'C': {
            'name': '2-styrylepoxide',
            'formula': 'C10H10O',  # C6H5-CH=CH-CH(O)CH2
            'aromatic_group': 'phenyl'
        },
        'D': {
            'name': '2-(4-methylstyryl)oxirane',
            'formula': 'C11H12O',  # CH3-C6H4-CH=CH-CH(O)CH2
            'aromatic_group': 'p-tolyl'
        }
    }

    llm_answer = 'D'

    # --- Step 3: Perform the verification ---
    # The reaction is an isomerization, so the starting material must have the same
    # molecular formula and the same core skeleton as the product.

    # Find all options that satisfy the constraints.
    valid_options = []
    for key, properties in options.items():
        # Constraint 1: Molecular formula must match C11H12O.
        formula_match = (properties['formula'] == required_formula)
        
        # Constraint 2: The aromatic group must match what's found in the product.
        group_match = (properties['aromatic_group'] == product_required_group)
        
        if formula_match and group_match:
            valid_options.append(key)

    # --- Step 4: Conclude based on the verification ---
    if len(valid_options) == 0:
        return "Incorrect. No option satisfies all constraints. The logic of the provided answer is flawed."
    
    if len(valid_options) > 1:
        return f"Incorrect. The problem is ambiguous as options {valid_options} all satisfy the constraints."

    # At this point, there is exactly one valid option. Check if it matches the LLM's answer.
    correct_option = valid_options[0]
    
    if llm_answer == correct_option:
        return "Correct"
    else:
        # This part explains why the LLM's answer would be wrong if it were different.
        llm_option_data = options[llm_answer]
        if llm_option_data['formula'] != required_formula:
            reason = f"its formula is {llm_option_data['formula']}, not {required_formula}."
        elif llm_option_data['aromatic_group'] != product_required_group:
            reason = f"it contains a {llm_option_data['aromatic_group']} group, but the product's NMR indicates a {product_required_group} group."
        else:
            reason = "of an unknown error." # Should not be reached
            
        return f"Incorrect. The provided answer '{llm_answer}' is wrong because {reason} The correct answer must be '{correct_option}'."

# Execute the check
result = check_chemistry_answer()
print(result)