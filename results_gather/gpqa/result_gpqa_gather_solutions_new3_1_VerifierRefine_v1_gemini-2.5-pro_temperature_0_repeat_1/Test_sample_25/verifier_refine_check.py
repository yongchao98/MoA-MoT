def check_chemistry_answer():
    """
    Checks the correctness of the selected answer for the given chemistry question.
    The question has two parts:
    A) Identifying a reactant in a [2+2] cycloaddition.
    B) Ordering the reactivity of dienes in a Diels-Alder reaction.
    """
    # Define the options from the question
    options = {
        'A': {
            'reactant_A': '2,2-diiodoethen-1-one',
            'reactivity_order_B': [3, 1, 2, 4]
        },
        'B': {
            'reactant_A': '4,4-diiodocyclobut-2-en-1-one',
            'reactivity_order_B': [4, 2, 1, 3]
        },
        'C': {
            'reactant_A': '4,4-diiodocyclobut-2-en-1-one',
            'reactivity_order_B': [3, 1, 2, 4]
        },
        'D': {
            'reactant_A': '2,2-diiodoethen-1-one',
            'reactivity_order_B': [4, 2, 1, 3]
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'A'

    # --- Correct answers based on chemical principles ---

    # Part A: Correct Reactant
    # The reaction is a [2+2] cycloaddition of cyclohexene with a ketene.
    # The product's substituents (ketone at C7, diiodo at C8) require the ketene to be I2C=C=O.
    # The IUPAC name for I2C=C=O is 2,2-diiodoethen-1-one.
    correct_reactant_A = '2,2-diiodoethen-1-one'

    # Part B: Correct Reactivity Order
    # Reactivity is based on the ability to form the s-cis conformation and electronic effects.
    # 3 (cyclopenta-1,3-diene): Locked s-cis, most reactive.
    # 1 (2,3-dimethylbuta-1,3-diene): Internal EDGs, very reactive.
    # 2 ((2E,4E)-hexa-2,4-diene): Terminal EDGs, reactive.
    # 4 ((2Z,4Z)-hexa-2,4-diene): Sterically hindered from forming s-cis, least reactive.
    # The correct order is 3 > 1 > 2 > 4.
    correct_reactivity_order_B = [3, 1, 2, 4]

    # --- Verification ---
    if llm_answer_choice not in options:
        return f"Error: The provided answer choice '{llm_answer_choice}' is not a valid option (A, B, C, or D)."

    selected_option = options[llm_answer_choice]

    # Check Part A
    if selected_option['reactant_A'] != correct_reactant_A:
        return (f"Incorrect. The reactant A in option {llm_answer_choice} is wrong.\n"
                f"Reason: The reaction is a [2+2] cycloaddition between cyclohexene and a ketene. "
                f"To form 8,8-diiodobicyclo[4.2.0]octan-7-one, the ketene must be diiodoketene (I2C=C=O), "
                f"which is named '{correct_reactant_A}'. The answer provided '{selected_option['reactant_A']}'.")

    # Check Part B
    if selected_option['reactivity_order_B'] != correct_reactivity_order_B:
        return (f"Incorrect. The reactivity order B in option {llm_answer_choice} is wrong.\n"
                f"Reason: The correct reactivity order (most to least) is {correct_reactivity_order_B}, based on "
                f"the diene's ability to adopt the s-cis conformation and electronic effects. The answer provided "
                f"{selected_option['reactivity_order_B']}.")

    # If both parts are correct
    return "Correct"

# Run the check and print the result
result = check_chemistry_answer()
print(result)