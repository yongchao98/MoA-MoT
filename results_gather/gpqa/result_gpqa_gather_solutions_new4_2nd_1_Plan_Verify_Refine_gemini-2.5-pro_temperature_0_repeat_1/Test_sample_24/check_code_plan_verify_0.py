def check_chemistry_answer():
    """
    Checks the correctness of the answer to a chemistry question about named reactions.
    """
    # The final answer provided by the LLM to be checked.
    llm_answer = 'D'

    # Define the four options with their proposed reactants.
    options = {
        'A': {
            'reactant_A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'reactant_B': '4-methyl-1-phenylpent-3-en-1-one'
        },
        'B': {
            'reactant_A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'reactant_B': '4-methyl-1-phenylpent-3-en-1-one'
        },
        'C': {
            'reactant_A': '2,8-dimethylspiro[4.5]decan-6-ol',
            'reactant_B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        },
        'D': {
            'reactant_A': '2,7-dimethyloctahydronaphthalene-4a,8a-diol',
            'reactant_B': '(((3-methylbut-2-en-1-yl)oxy)methyl)benzene'
        }
    }

    # --- Define the logical checks based on chemical principles ---

    # Rule for Reaction A (Pinacol Rearrangement):
    # The reaction converts a 1,2-diol into a ketone.
    # The reactant must be a diol. We can check this by looking for "diol" in the name.
    def is_valid_for_pinacol(reactant_name):
        return "diol" in reactant_name

    # Rule for Reaction B (Wittig Rearrangement):
    # The reaction converts an ether into a rearranged alcohol using a strong base.
    # The reactant must be an ether. We can identify ethers by the "oxy" in their name.
    # It cannot be a ketone ("-one") because BuLi would add a butyl group, which is not in the product.
    def is_valid_for_wittig(reactant_name):
        return "oxy" in reactant_name and "one" not in reactant_name

    # --- Apply the checks to the provided answer ---
    
    # Retrieve the reactants from the LLM's chosen option
    chosen_option_reactants = options.get(llm_answer)
    
    if not chosen_option_reactants:
        return f"Incorrect. The provided answer '{llm_answer}' is not one of the valid options (A, B, C, D)."

    # Check Reactant A from the chosen option
    if not is_valid_for_pinacol(chosen_option_reactants['reactant_A']):
        return (f"Incorrect. The provided answer is {llm_answer}. "
                f"Its reactant A, '{chosen_option_reactants['reactant_A']}', is not a diol. "
                f"A Pinacol rearrangement requires a 1,2-diol as the starting material.")

    # Check Reactant B from the chosen option
    if not is_valid_for_wittig(chosen_option_reactants['reactant_B']):
        return (f"Incorrect. The provided answer is {llm_answer}. "
                f"Its reactant B, '{chosen_option_reactants['reactant_B']}', is not an ether. "
                f"A Wittig rearrangement requires an ether as the starting material.")

    # --- Verify that the chosen answer is the *only* correct one ---
    
    correctly_identified_options = []
    for option_key, reactants in options.items():
        # Check if both reactants in the current option satisfy the rules
        if is_valid_for_pinacol(reactants['reactant_A']) and is_valid_for_wittig(reactants['reactant_B']):
            correctly_identified_options.append(option_key)
            
    if len(correctly_identified_options) > 1:
        return (f"Incorrect. The logic leads to multiple correct options: {correctly_identified_options}. "
                f"The question may be ambiguous.")
    
    if not correctly_identified_options:
        return "Incorrect. The logic shows that none of the options are fully correct."

    # Compare the LLM's answer with the logically derived correct answer
    if llm_answer == correctly_identified_options[0]:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but a systematic check reveals "
                f"the correct answer should be {correctly_identified_options[0]}.")

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)