def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer by applying chemical principles.
    1. Determines the correct stereochemistry for starting material A based on LiBH4 selectivity and CIP priority changes.
    2. Determines the correct stereochemistry for starting material B based on BH3 selectivity and CIP priority changes.
    3. Compares the derived correct answer with the LLM's proposed answer.
    """
    # The final answer provided by the LLM to be checked.
    llm_answer_choice = 'B'

    # Define the options as per the question.
    options = {
        'A': {'A': '(R)', 'B': '(S)'},
        'B': {'A': '(S)', 'B': '(S)'},
        'C': {'A': '(S)', 'B': '(R)'},
        'D': {'A': '(R)', 'B': '(R)'}
    }

    # --- Step 1: Logic for Reaction A ---
    # Reaction: A + LiBH4 -> (R)-product
    # LiBH4 reduces the ester group (-CH2COOiBu) to an alcohol (-CH2CH2OH).
    # The CIP priority of the side chains on the chiral center changes from:
    # (-CH2COOiBu > -CH2COOH) to (-CH2COOH > -CH2CH2OH).
    # This swap of the top two priorities inverts the R/S designation.
    # Therefore, to get an (R) product, the starting material A must be (S).
    correct_A_config = '(S)'

    # --- Step 2: Logic for Reaction B ---
    # Reaction: B + BH3 -> (S)-product
    # BH3 reduces the carboxylic acid group (-CH2COOH) to an alcohol (-CH2CH2OH).
    # The CIP priority of the side chains on the chiral center changes from:
    # (-CH2COOiBu > -CH2COOH) to (-CH2COOiBu > -CH2CH2OH).
    # The priority order does NOT change. This retains the R/S designation.
    # Therefore, to get an (S) product, the starting material B must be (S).
    correct_B_config = '(S)'

    # --- Step 3: Check the LLM's answer ---
    if llm_answer_choice not in options:
        return f"Invalid Answer Format: The answer '{llm_answer_choice}' is not a valid option (A, B, C, or D)."

    proposed_config = options[llm_answer_choice]
    errors = []

    # Check constraint for A
    if proposed_config['A'] != correct_A_config:
        errors.append(
            f"Constraint for starting material A is not satisfied. "
            f"The reaction A + LiBH4 -> (R)-product requires A to be {correct_A_config}, but the answer proposes {proposed_config['A']}. "
            f"Reason: The reduction of the ester with LiBH4 causes an inversion of the R/S designation due to a change in CIP priorities."
        )

    # Check constraint for B
    if proposed_config['B'] != correct_B_config:
        errors.append(
            f"Constraint for starting material B is not satisfied. "
            f"The reaction B + BH3 -> (S)-product requires B to be {correct_B_config}, but the answer proposes {proposed_config['B']}. "
            f"Reason: The reduction of the carboxylic acid with BH3 does not change the CIP priority order, resulting in retention of the R/S designation."
        )

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Execute the check
result = check_chemistry_answer()
print(result)