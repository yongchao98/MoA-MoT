def check_answer():
    """
    Checks the correctness of the final answer based on chemical principles.
    """
    # Define the options as provided in the question
    options = {
        "A": {
            "reactant_A": "4,4-diiodocyclobut-2-en-1-one",
            "reactivity_B": [3, 1, 2, 4]
        },
        "B": {
            "reactant_A": "4,4-diiodocyclobut-2-en-1-one",
            "reactivity_B": [4, 2, 1, 3]
        },
        "C": {
            "reactant_A": "2,2-diiodoethen-1-one",
            "reactivity_B": [4, 2, 1, 3]
        },
        "D": {
            "reactant_A": "2,2-diiodoethen-1-one",
            "reactivity_B": [3, 1, 2, 4]
        }
    }

    # The final answer to be checked
    final_answer_key = "D"

    # Define the correct answer based on chemical analysis
    correct_reactant_A = "2,2-diiodoethen-1-one"
    correct_reactivity_B = [3, 1, 2, 4]

    # Retrieve the selected answer's content
    selected_option = options.get(final_answer_key)

    if not selected_option:
        return f"Invalid answer key '{final_answer_key}'. The key must be one of {list(options.keys())}."

    # Check if the selected option's parts match the correct analysis
    is_A_correct = selected_option["reactant_A"] == correct_reactant_A
    is_B_correct = selected_option["reactivity_B"] == correct_reactivity_B

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            reason_A = (
                "Constraint on reactant A is not satisfied. The reaction is a [2+2] cycloaddition "
                "requiring diiodoketene (I₂C=C=O), which is named '2,2-diiodoethen-1-one'. "
                f"The answer incorrectly identifies reactant A as '{selected_option['reactant_A']}'."
            )
            error_messages.append(reason_A)
        
        if not is_B_correct:
            reason_B = (
                "Constraint on diene reactivity order B is not satisfied. The reactivity (most>least) is determined by "
                "the ability to form the s-cis conformation. The correct order is 3 (locked s-cis) > 1 (easy s-cis, EDGs) > "
                "2 (possible s-cis) > 4 (hindered s-cis). "
                f"The correct sequence is {correct_reactivity_B}, but the answer chose {selected_option['reactivity_B']}."
            )
            error_messages.append(reason_B)
        
        return "\n".join(error_messages)

# Run the check and print the result
print(check_answer())