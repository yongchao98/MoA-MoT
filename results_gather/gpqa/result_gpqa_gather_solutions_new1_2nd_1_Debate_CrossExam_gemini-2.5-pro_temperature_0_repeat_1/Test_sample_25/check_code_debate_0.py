def check_organic_chemistry_answer():
    """
    This function checks the correctness of the selected answer for the given organic chemistry question.
    It verifies two parts:
    1. The identity of reactant A.
    2. The reactivity order of the dienes in B.
    """

    # Define the options as provided in the question prompt.
    # The keys 'A' and 'B' in the inner dictionaries correspond to the two parts of the answer.
    options = {
        'A': {
            'A': "4,4-diiodocyclobut-2-en-1-one",
            'B': [4, 2, 1, 3]
        },
        'B': {
            'A': "4,4-diiodocyclobut-2-en-1-one",
            'B': [3, 1, 2, 4]
        },
        'C': {
            'A': "2,2-diiodoethen-1-one",
            'B': [4, 2, 1, 3]
        },
        'D': {
            'A': "2,2-diiodoethen-1-one",
            'B': [3, 1, 2, 4]
        }
    }

    # The final answer to be checked, as provided by the LLM.
    llm_final_answer_key = "D"

    # --- Part 1: Determine the correct identity of Reactant A ---
    # The reaction is: Cyclohexene + A ---> 8,8-diiodobicyclo[4.2.0]octan-7-one.
    # This is a [2+2] cycloaddition. The product is a fusion of a 6-membered ring (from cyclohexene)
    # and a 4-membered ring (from A).
    # The 4-membered ring contains a ketone (C=O) and a diiodo-substituted carbon (CI2).
    # This requires reactant A to be a ketene with the structure I2C=C=O.
    # The systematic IUPAC name for I2C=C=O is "2,2-diiodoethen-1-one".
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # --- Part 2: Determine the correct reactivity order for B ---
    # The reactivity of dienes in cycloaddition reactions depends on:
    # 1. Ability to adopt s-cis conformation: Dienes locked in s-cis are most reactive.
    #    Dienes sterically hindered from forming s-cis are least reactive.
    # 2. Electronic effects: Electron-donating groups (EDGs) increase reactivity.
    #
    # Analysis of the dienes:
    # 3. cyclopenta-1,3-diene: Locked in s-cis conformation. MOST REACTIVE.
    # 4. (2Z,4Z)-hexa-2,4-diene: Severe steric hindrance prevents s-cis conformation. LEAST REACTIVE.
    # 1. 2,3-dimethylbuta-1,3-diene: Two internal EDGs, highly activating. More reactive than 2.
    # 2. (2E,4E)-hexa-2,4-diene: Two terminal EDGs, less activating than internal ones. Less reactive than 1.
    #
    # The correct order from most to least reactive is 3 > 1 > 2 > 4.
    correct_reactivity_order_B = [3, 1, 2, 4]

    # --- Part 3: Validate the selected answer ---
    if llm_final_answer_key not in options:
        return f"Error: The provided answer key '{llm_final_answer_key}' is not a valid option (A, B, C, or D)."

    selected_option = options[llm_final_answer_key]
    
    # Check if both parts of the selected answer match the correct analysis
    is_A_correct = (selected_option['A'] == correct_reactant_A)
    is_B_correct = (selected_option['B'] == correct_reactivity_order_B)

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            error_messages.append(
                f"Constraint on Reactant A is not satisfied. "
                f"The answer claims A is '{selected_option['A']}', but the correct reactant is '{correct_reactant_A}'."
            )
        if not is_B_correct:
            error_messages.append(
                f"Constraint on Reactivity Order B is not satisfied. "
                f"The answer claims the order is {selected_option['B']}, but the correct order is {correct_reactivity_order_B}."
            )
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_organic_chemistry_answer()
print(result)