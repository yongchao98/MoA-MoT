def check_answer():
    """
    Checks the correctness of the LLM's answer for the given chemistry question.
    """
    # Correct product names based on chemical principles
    correct_product_A = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"
    correct_product_B = "3-methyl-4-nitrohexanenitrile"

    # Options provided in the question
    options = {
        "A": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "B": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        },
        "C": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "D": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        }
    }

    # The final answer provided by the LLM
    llm_answer_choice = "B"

    # Retrieve the product names from the selected option
    selected_option = options.get(llm_answer_choice)

    if not selected_option:
        return f"Invalid answer choice '{llm_answer_choice}'. Please choose from A, B, C, or D."

    # Check if the names in the selected option match the correct names
    is_A_correct = selected_option["A"] == correct_product_A
    is_B_correct = selected_option["B"] == correct_product_B

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            error_messages.append(
                f"Product A is incorrect. The selected answer provides '{selected_option['A']}', but the correct product is '{correct_product_A}'. The error is in the naming of the cyclohexyl substituent, specifically the position of the oxo group and methyl groups."
            )
        if not is_B_correct:
            error_messages.append(
                f"Product B is incorrect. The selected answer provides '{selected_option['B']}', but the correct product is '{correct_product_B}'. The error is in the parent chain length and substituent positions."
            )
        return "\n".join(error_messages)

# Run the check
result = check_answer()
print(result)