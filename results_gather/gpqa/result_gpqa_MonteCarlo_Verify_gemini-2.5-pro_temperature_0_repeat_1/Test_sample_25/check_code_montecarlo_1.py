def check_chemistry_answer():
    """
    This function programmatically checks the correctness of the selected answer
    for the given organic chemistry question. It validates both the identification
    of reactant A and the reactivity order of the dienes in B.
    """

    # --- Ground Truth Analysis ---

    # Part A: Identification of Reactant A
    # The reaction is a [2+2] cycloaddition between cyclohexene and a ketene.
    # The product, 8,8-diiodobicyclo[4.2.0]octan-7-one, dictates that the
    # ketene must be I2C=C=O.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # Part B: Diene Reactivity Order (for Diels-Alder reactions)
    # 1. 2,3-dimethylbuta-1,3-diene: Acyclic, internal EDGs (activating).
    # 2. (2E,4E)-hexa-2,4-diene: Acyclic, terminal EDGs (less activating than internal).
    # 3. cyclopenta-1,3-diene: Cyclic, locked in reactive s-cis conformation (most reactive).
    # 4. (2Z,4Z)-hexa-2,4-diene: Acyclic, severe steric hindrance prevents s-cis (least reactive).
    # The order from most to least reactive is 3 > 1 > 2 > 4.
    correct_reactivity_order_B = [3, 1, 2, 4]

    # --- Options from the Question ---
    options = {
        'A': {'A': '4,4-diiodocyclobut-2-en-1-one', 'B': [3, 1, 2, 4]},
        'B': {'A': '4,4-diiodocyclobut-2-en-1-one', 'B': [4, 2, 1, 3]},
        'C': {'A': '2,2-diiodoethen-1-one', 'B': [3, 1, 2, 4]},
        'D': {'A': '2,2-diiodoethen-1-one', 'B': [4, 2, 1, 3]}
    }

    # The answer to be checked is 'C'
    answer_to_check = 'C'
    selected_option = options[answer_to_check]

    # --- Verification Logic ---
    is_A_correct = (selected_option['A'] == correct_reactant_A)
    is_B_correct = (selected_option['B'] == correct_reactivity_order_B)

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            error_messages.append(
                f"Constraint for reactant A is not satisfied. "
                f"Expected: '{correct_reactant_A}', but got: '{selected_option['A']}'. "
                f"Reason: The reaction is a [2+2] cycloaddition forming a bicyclo[4.2.0] system with a ketone and diiodo group, which requires diiodoketene as the reactant."
            )
        if not is_B_correct:
            error_messages.append(
                f"Constraint for diene reactivity order B is not satisfied. "
                f"Expected: {correct_reactivity_order_B}, but got: {selected_option['B']}. "
                f"Reason: The correct order is based on s-cis conformation availability (cyclopentadiene is locked s-cis) and electronic effects (internal EDGs are more activating)."
            )
        return "\n".join(error_messages)

# Execute the checker and print the result
result = check_chemistry_answer()
print(result)