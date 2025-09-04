def check_chemistry_answer():
    """
    This function checks the correctness of the selected answer based on two key chemical principles
    of the Michael addition reactions described in the question.
    """

    # The options provided in the question
    options = {
        'A': {
            'A': "methyl 3-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            'B': "4-ethyl 1-methyl 2-cyclopentyl-3,3-diethyl-2-phenylsuccinate"
        },
        'B': {
            'A': "methyl 1-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            'B': "4-ethyl 1-methyl 2-cyclopentyl-3,3-diethyl-2-phenylsuccinate"
        },
        'C': {
            'A': "methyl 3-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            'B': "ethyl 2-ethyl-2-(1-(2-methoxy-2-oxo-1-phenylethyl)cyclopentyl)butanoate"
        },
        'D': {
            'A': "methyl 1-(2-((2,4-dimethylphenyl)sulfinyl)ethyl)-2-oxocyclohexane-1-carboxylate",
            'B': "ethyl 2-ethyl-2-(1-(2-methoxy-2-oxo-1-phenylethyl)cyclopentyl)butanoate"
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer = 'D'

    # --- Constraint 1: Check the regioselectivity of Reaction A ---
    # The most acidic proton is at C1 (between the two carbonyls).
    # Therefore, the product must be substituted at the C1 position.
    # The correct name must contain "methyl 1-".
    def check_product_A(name):
        return "methyl 1-" in name

    # --- Constraint 2: Check the core structure of Product B ---
    # The product is a Michael adduct of a butanoate enolate.
    # The product must be a substituted butanoate, not a succinate.
    def check_product_B(name):
        return "butanoate" in name and "succinate" not in name

    # Get the product names for the selected answer
    answer_to_check = options.get(llm_answer)

    if not answer_to_check:
        return f"Invalid answer key: '{llm_answer}'. The key must be one of A, B, C, or D."

    # Verify the selected answer against the chemical principles
    is_A_correct = check_product_A(answer_to_check['A'])
    if not is_A_correct:
        return (f"The provided answer '{llm_answer}' is incorrect. "
                f"Constraint for Reaction A is not satisfied. "
                f"Reason: The Michael addition on methyl 2-oxocyclohexane-1-carboxylate should occur at the C1 position, "
                f"resulting in a 'methyl 1-(...)' product. The answer's product A is named '{answer_to_check['A']}', which indicates incorrect regioselectivity.")

    is_B_correct = check_product_B(answer_to_check['B'])
    if not is_B_correct:
        return (f"The provided answer '{llm_answer}' is incorrect. "
                f"Constraint for Reaction B is not satisfied. "
                f"Reason: The Michael addition product should be a substituted butanoate, not a succinate. "
                f"The answer's product B is named '{answer_to_check['B']}', which is a structurally incorrect class of compound for this reaction.")

    # Final check to ensure no other option is also correct
    correct_options_found = []
    for option_key, products in options.items():
        if check_product_A(products['A']) and check_product_B(products['B']):
            correct_options_found.append(option_key)
    
    if len(correct_options_found) == 1 and correct_options_found[0] == llm_answer:
        return "Correct"
    else:
        return (f"The provided answer '{llm_answer}' is incorrect. "
                f"Based on the chemical principles, the only correct option is '{correct_options_found[0]}'.")

# Execute the check and print the result
result = check_chemistry_answer()
print(result)