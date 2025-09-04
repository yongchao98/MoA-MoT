def check_answer():
    """
    This function checks the correctness of the final answer for the given chemistry question.
    It verifies both parts of the answer: the identity of reactant A and the reactivity order of dienes B.
    """

    # --- Part 1: Define the correct answer based on chemical principles ---

    # Correct identity of reactant A:
    # The reaction is a [2+2] cycloaddition between cyclohexene and a ketene.
    # The product, 8,8-diiodobicyclo[4.2.0]octan-7-one, requires the ketene to be
    # diiodoketene (I2C=C=O), which is named 2,2-diiodoethen-1-one.
    correct_A = "2,2-diiodoethen-1-one"

    # Correct reactivity order for B (most reactive to least reactive):
    # 3. cyclopenta-1,3-diene: Locked in s-cis conformation, most reactive.
    # 1. 2,3-dimethylbuta-1,3-diene: Electron-donating groups (EDGs) on internal carbons, easily forms s-cis.
    # 2. (2E,4E)-hexa-2,4-diene: EDGs on terminal carbons, can form s-cis but less favored than (1).
    # 4. (2Z,4Z)-hexa-2,4-diene: Severe steric hindrance prevents s-cis conformation, least reactive.
    # The correct order is 3 > 1 > 2 > 4.
    correct_B_order = [3, 1, 2, 4]

    # --- Part 2: Define the options given in the question ---
    options = {
        "A": {
            "A": "2,2-diiodoethen-1-one",
            "B": [4, 2, 1, 3]
        },
        "B": {
            "A": "4,4-diiodocyclobut-2-en-1-one",
            "B": [4, 2, 1, 3]
        },
        "C": {
            "A": "4,4-diiodocyclobut-2-en-1-one",
            "B": [3, 1, 2, 4]
        },
        "D": {
            "A": "2,2-diiodoethen-1-one",
            "B": [3, 1, 2, 4]
        }
    }

    # --- Part 3: Extract and check the proposed answer ---
    
    # The proposed answer is 'D'.
    proposed_answer_key = "D"
    
    if proposed_answer_key not in options:
        return f"Invalid answer key '{proposed_answer_key}'. The valid options are A, B, C, D."

    proposed_option = options[proposed_answer_key]

    # Check if reactant A in the proposed option is correct.
    if proposed_option["A"] != correct_A:
        return (f"Incorrect. The proposed answer '{proposed_answer_key}' is wrong because the identity of reactant A is incorrect. "
                f"The option states A is '{proposed_option['A']}', but the correct reactant is '{correct_A}'.")

    # Check if the reactivity order B in the proposed option is correct.
    if proposed_option["B"] != correct_B_order:
        return (f"Incorrect. The proposed answer '{proposed_answer_key}' is wrong because the reactivity order B is incorrect. "
                f"The option states the order is {proposed_option['B']}, but the correct order is {correct_B_order}.")

    # If both parts are correct, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_answer()
print(result)