def check_answer():
    """
    Checks the correctness of the LLM's answer regarding optical isomerism.
    """
    # Ground truth based on chemical principles.
    # A compound shows optical isomerism if it is chiral.
    compound_properties = {
        1: {
            "name": "dimethyl 6,6'-dinitro-[1,1'-biphenyl]-2,2'-dicarboxylate",
            "is_chiral": True,
            "reason": "Exhibits atropisomerism due to restricted rotation around the biphenyl bond, making it chiral."
        },
        2: {
            "name": "methyl 2-hydroxypropanoate",
            "is_chiral": True,
            "reason": "Contains a chiral carbon center bonded to four different groups (-H, -OH, -CH3, -COOCH3)."
        },
        3: {
            "name": "benzophenone",
            "is_chiral": False,
            "reason": "Is achiral because it has a plane of symmetry."
        },
        4: {
            "name": "dimethyl fumarate",
            "is_chiral": False,
            "reason": "Is achiral because it is planar and has a center of inversion."
        }
    }

    # Define the options from the question
    options = {
        "A": {3, 4},
        "B": {2, 3},
        "C": {1, 2, 4},
        "D": {1, 2}
    }

    # The final answer provided by the LLM
    llm_answer_choice = "D"

    # Determine the correct set of optically active compounds
    correct_chiral_set = {num for num, prop in compound_properties.items() if prop["is_chiral"]}

    # Find the correct option letter corresponding to the correct set
    correct_option = None
    for letter, compound_set in options.items():
        if compound_set == correct_chiral_set:
            correct_option = letter
            break

    # Check if the LLM's answer matches the correct option
    if llm_answer_choice == correct_option:
        return "Correct"
    else:
        error_message = f"Incorrect. The provided answer is {llm_answer_choice}, but the correct answer is {correct_option}.\n"
        error_message += "Reasoning:\n"
        for i in range(1, 5):
            is_chiral_text = "shows" if compound_properties[i]["is_chiral"] else "does not show"
            error_message += f"- Compound {i} ({compound_properties[i]['name']}) {is_chiral_text} optical isomerism. Reason: {compound_properties[i]['reason']}\n"
        error_message += f"The set of compounds showing optical isomerism is {correct_chiral_set}, which corresponds to option {correct_option}."
        return error_message

# Run the check
result = check_answer()
print(result)