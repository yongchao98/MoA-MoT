def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by applying chemical principles.
    """

    # --- Define Chemical Principles ---

    # 1. Reagent Selectivity: Which functional group is reduced?
    # The starting material has a carboxylic acid and an ester.
    # LiBH4 reduces the ester. BH3 reduces the carboxylic acid.
    reagent_selectivity = {
        "LiBH4": "ester",
        "BH3": "carboxylic_acid"
    }

    # 2. Stereochemistry Principle
    # The reactions do not involve the chiral center. Thus, the configuration is retained.
    stereochemistry_rule = "retained"

    # --- Define the Problem Statement ---

    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    reactions = {
        "A": {"reagent": "LiBH4", "product_config": "R"},
        "B": {"reagent": "BH3", "product_config": "S"}
    }

    # --- Define the LLM's Answer ---
    llm_answer_option = 'B'
    options = {
        'A': {'A': 'S', 'B': 'R'},
        'B': {'A': 'R', 'B': 'S'},
        'C': {'A': 'S', 'B': 'S'},
        'D': {'A': 'R', 'B': 'R'}
    }
    llm_proposed_configs = options.get(llm_answer_option)

    if llm_proposed_configs is None:
        return f"The provided answer '{llm_answer_option}' is not a valid option."

    # --- Verification Process ---

    # Step 1: Deduce the required configuration for each starting material.
    deduced_configs = {}
    error_messages = []

    for material, reaction_info in reactions.items():
        product_config = reaction_info["product_config"]
        
        # Apply the stereochemistry rule
        if stereochemistry_rule == "retained":
            # If the configuration is retained, the starting material must have the same configuration as the product.
            required_start_config = product_config
        elif stereochemistry_rule == "inverted":
            required_start_config = "S" if product_config == "R" else "R"
        else:
            # This case handles potential errors or other unstated rules like racemization.
            return f"Error: The stereochemistry rule '{stereochemistry_rule}' is not recognized."
        
        deduced_configs[material] = required_start_config

    # Step 2: Compare the deduced configurations with the LLM's answer.
    if deduced_configs == llm_proposed_configs:
        return "Correct"
    else:
        # Generate a detailed error report if the answer is incorrect.
        reason = "The provided answer is incorrect because the deduced starting material configurations do not match the chosen option.\n"
        reason += f"Based on the principle that stereochemistry is retained:\n"
        
        # Check material A
        if deduced_configs['A'] != llm_proposed_configs['A']:
            reason += f"- For Reaction A, to produce the (R)-product, the starting material A must have the (R) configuration. The answer incorrectly proposes it is ({llm_proposed_configs['A']}).\n"
        
        # Check material B
        if deduced_configs['B'] != llm_proposed_configs['B']:
            reason += f"- For Reaction B, to produce the (S)-product, the starting material B must have the (S) configuration. The answer incorrectly proposes it is ({llm_proposed_configs['B']}).\n"
            
        reason += f"\nCorrect Deduction: A should be ({deduced_configs['A']}) and B should be ({deduced_configs['B']})."
        reason += f"\nProposed Answer (Option {llm_answer_option}): A is ({llm_proposed_configs['A']}) and B is ({llm_proposed_configs['B']})."
        return reason

# Run the check
result = check_chemistry_answer()
print(result)