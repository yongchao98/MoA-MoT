def check_chemistry_answer():
    """
    Checks the correctness of the selected answer for the given organic chemistry question.
    """
    # Based on chemical principles, the correct names for the products are determined.
    correct_product_A_name = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"
    correct_product_B_name = "3-methyl-4-nitrohexanenitrile"

    # The options provided in the question.
    options = {
        'A': {
            'A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'B': "3-methyl-4-nitrohexanenitrile"
        },
        'B': {
            'A': "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            'B': "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'C': {
            'A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'B': "2,3-dimethyl-4-nitrobutanenitrile"
        },
        'D': {
            'A': "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            'B': "3-methyl-4-nitrohexanenitrile"
        }
    }

    # The final answer provided by the LLM.
    llm_answer_key = "A"

    # Get the proposed names from the selected option.
    selected_option = options.get(llm_answer_key)
    if not selected_option:
        return f"Error: The answer key '{llm_answer_key}' is not a valid option (A, B, C, or D)."

    proposed_A = selected_option['A']
    proposed_B = selected_option['B']

    # Check if the proposed names match the correct names.
    errors = []
    if proposed_A != correct_product_A_name:
        reason = (f"Product A is incorrect. The proposed name is '{proposed_A}'.\n"
                  f"The correct name should be '{correct_product_A_name}'.\n"
                  f"Reason: The Michael addition occurs at the C6 position of the ketone (the only position with an alpha-proton). "
                  f"When naming the product as an ethyl propanoate derivative, the cyclohexyl substituent is numbered from the point of attachment, "
                  f"leading to the name (3-ethyl-1,3-dimethyl-2-oxocyclohexyl).")
        errors.append(reason)

    if proposed_B != correct_product_B_name:
        reason = (f"Product B is incorrect. The proposed name is '{proposed_B}'.\n"
                  f"The correct name should be '{correct_product_B_name}'.\n"
                  f"Reason: The Michael addition between 1-nitropropane and but-2-enenitrile results in a product with a 6-carbon main chain (hexanenitrile). "
                  f"Numbering from the nitrile group places the methyl group at C3 and the nitro group at C4.")
        errors.append(reason)

    if not errors:
        return "Correct"
    else:
        return "Incorrect. " + "\n\n".join(errors)

# Execute the check
result = check_chemistry_answer()
print(result)