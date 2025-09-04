import re

def check_answer():
    """
    Checks the correctness of the final answer based on the analysis of the two chemical reactions.
    """
    # Data from the question options
    options = {
        'A': {
            'product_A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'product_B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'B': {
            'product_A': "4-methyl-1-phenylpent-3-en-1-ol",
            'product_B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        },
        'C': {
            'product_A': "4-methyl-1-phenylpent-3-en-1-ol",
            'product_B': "2,3,4,6,7,8-hexamethyl-5,9,10,11-tetrahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorine"
        },
        'D': {
            'product_A': "(Z)-2-methyl-5-phenylpent-2-en-1-ol",
            'product_B': "2,3,4,6,7,8-hexamethyl-5,9,9a,10,11,11a-hexahydro-1H-benzo[3,4]azuleno[1,8,7,6-cdef]fluorene"
        }
    }

    # The final answer provided by the LLM
    llm_answer_text = "<<<B>>>"
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."
    
    llm_choice = match.group(1)

    # --- Define Constraints based on Chemical Principles ---

    # Constraint 1: Analysis of Reaction A (Wittig Rearrangement)
    # The most plausible product from the [1,2]-Wittig rearrangement is 4-methyl-1-phenylpent-3-en-1-ol.
    # The other option, (Z)-2-methyl-5-phenylpent-2-en-1-ol, is a primary alcohol with a different skeleton, which is not a direct product of this rearrangement.
    expected_product_A = "4-methyl-1-phenylpent-3-en-1-ol"

    # Constraint 2: Analysis of Reaction B (Cope Rearrangement)
    # The Cope rearrangement is an isomerization. The product must have the same molecular formula and degree of saturation as the starting material.
    # The starting material is a "hexahydro" derivative, so the product must also be a "hexahydro" derivative.
    expected_saturation_B = "hexahydro"

    # --- Check all options against constraints to find the correct one(s) ---
    correct_options = []
    for option, products in options.items():
        # Check Constraint 1
        is_A_correct = (products['product_A'] == expected_product_A)
        
        # Check Constraint 2
        is_B_correct = (expected_saturation_B in products['product_B'])
        
        if is_A_correct and is_B_correct:
            correct_options.append(option)

    # --- Validate the LLM's answer ---
    if not correct_options:
        return "Logic Error: None of the options satisfy both chemical constraints. There might be an error in the question or the constraints."
    
    if len(correct_options) > 1:
        return f"Logic Error: The constraints are not sufficient to uniquely determine the answer. Options {correct_options} are all valid."

    # The single correct option based on the logic
    true_correct_option = correct_options[0]

    if llm_choice == true_correct_option:
        return "Correct"
    else:
        # Provide a specific reason for the failure
        chosen_option_data = options[llm_choice]
        
        # Check why the chosen option is wrong
        is_A_correct_for_choice = (chosen_option_data['product_A'] == expected_product_A)
        is_B_correct_for_choice = (expected_saturation_B in chosen_option_data['product_B'])

        reasons = []
        if not is_A_correct_for_choice:
            reasons.append(f"Constraint 1 (Reaction A) is not satisfied. The product A in option {llm_choice} is '{chosen_option_data['product_A']}', but the expected product is '{expected_product_A}'.")
        
        if not is_B_correct_for_choice:
            reasons.append(f"Constraint 2 (Reaction B) is not satisfied. The product B in option {llm_choice} is a '{'tetrahydro' if 'tetrahydro' in chosen_option_data['product_B'] else 'unknown saturation'}' derivative, but the reaction is an isomerization and should produce a '{expected_saturation_B}' derivative.")
            
        return f"Incorrect: The chosen answer '{llm_choice}' is wrong. The correct answer should be '{true_correct_option}'.\nReason(s): {' '.join(reasons)}"

# Run the check
result = check_answer()
print(result)