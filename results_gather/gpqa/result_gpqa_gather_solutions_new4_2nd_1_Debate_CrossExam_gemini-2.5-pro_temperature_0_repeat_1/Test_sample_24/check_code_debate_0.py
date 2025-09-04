def check_chemistry_answer():
    """
    Checks the correctness of the answer to the chemistry question.
    The function encodes the chemical principles for each reaction and verifies
    if the provided answer's reactants satisfy these principles.
    """

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    # Define the options from the question.
    options = {
        "A": {
            "reactant_A": "2,8-dimethylspiro[4.5]decan-6-ol",
            "reactant_B": "4-methyl-1-phenylpent-3-en-1-one"
        },
        "B": {
            "reactant_A": "2,7-dimethyloctahydronaphthalene-4a,8a-diol",
            "reactant_B": "(((3-methylbut-2-en-1-yl)oxy)methyl)benzene"
        },
        "C": {
            "reactant_A": "2,8-dimethylspiro[4.5]decan-6-ol",
            "reactant_B": "(((3-methylbut-2-en-1-yl)oxy)methyl)benzene"
        },
        "D": {
            "reactant_A": "2,7-dimethyloctahydronaphthalene-4a,8a-diol",
            "reactant_B": "4-methyl-1-phenylpent-3-en-1-one"
        }
    }

    # --- Verification Logic ---

    # Rule for Reaction A (Pinacol Rearrangement): Reactant must be a 1,2-diol.
    # We can check this by seeing if the name contains "-diol".
    def check_reactant_A(name):
        return "diol" in name

    # Rule for Reaction B (Wittig Rearrangement): Reactant must be an ether.
    # We can check this by seeing if the name contains "-oxy-" and is not a ketone ("-one").
    def check_reactant_B(name):
        return "oxy" in name and "one" not in name

    # Determine the theoretically correct option based on the rules.
    correct_option = None
    for option_key, reactants in options.items():
        if check_reactant_A(reactants["reactant_A"]) and check_reactant_B(reactants["reactant_B"]):
            correct_option = option_key
            break
    
    # Compare the LLM's answer with the derived correct option.
    if llm_answer == correct_option:
        return "Correct"
    else:
        # If the answer is incorrect, provide a specific reason.
        chosen_reactants = options.get(llm_answer)
        if not chosen_reactants:
            return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

        # Check which constraint is violated by the given answer.
        if not check_reactant_A(chosen_reactants["reactant_A"]):
            return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option}.\n"
                    f"Reason: Reaction A is a Pinacol Rearrangement, which requires a 1,2-diol as the starting material. "
                    f"Option {llm_answer} incorrectly proposes '{chosen_reactants['reactant_A']}', which is an alcohol, not a diol.")
        
        if not check_reactant_B(chosen_reactants["reactant_B"]):
            return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option}.\n"
                    f"Reason: Reaction B is a Wittig Rearrangement, which requires an ether as the starting material. "
                    f"Option {llm_answer} incorrectly proposes '{chosen_reactants['reactant_B']}', which is a ketone, not an ether.")
        
        return f"Incorrect. The provided answer {llm_answer} is wrong for an undetermined reason, but the correct answer is {correct_option}."

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)