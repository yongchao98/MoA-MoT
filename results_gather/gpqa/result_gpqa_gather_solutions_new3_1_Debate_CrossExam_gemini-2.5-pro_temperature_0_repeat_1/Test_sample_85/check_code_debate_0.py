def check_chemistry_answer():
    """
    Checks the correctness of the answer to the chemistry question.

    The logic is based on two principles:
    1. Stereochemistry Retention: The reactions (reduction and cyclization) do not affect the chiral center.
       Therefore, an (R)-starting material yields an (R)-product, and an (S)-starting material yields an (S)-product.
    2. Reagent Selectivity:
       - LiBH4 selectively reduces the ester group.
       - BH3 selectively reduces the carboxylic acid group.
       Both reaction pathways lead to an intermediate that can cyclize to the final lactone product.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = "D"

    # Define the problem from the question
    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    product_A_stereochem = 'R'
    product_B_stereochem = 'S'

    # --- Step 1: Determine the required stereochemistry for starting material A ---
    # Principle: Stereochemistry is retained.
    # To get an (R)-product from Reaction A, starting material A must be (R).
    required_A_stereochem = product_A_stereochem

    # --- Step 2: Determine the required stereochemistry for starting material B ---
    # Principle: Stereochemistry is retained.
    # To get an (S)-product from Reaction B, starting material B must be (S).
    required_B_stereochem = product_B_stereochem

    # --- Step 3: Define the options provided in the question ---
    options = {
        "A": {"A": "R", "B": "R"},
        "B": {"A": "S", "B": "S"},
        "C": {"A": "S", "B": "R"},
        "D": {"A": "R", "B": "S"}
    }

    # --- Step 4: Check if the LLM's answer matches the derived correct configuration ---
    if llm_final_answer not in options:
        return f"Invalid answer format. The answer '{llm_final_answer}' is not one of the options A, B, C, or D."

    chosen_option_config = options[llm_final_answer]

    # Check if the configuration for A in the chosen option is correct
    is_A_correct = (chosen_option_config["A"] == required_A_stereochem)
    # Check if the configuration for B in the chosen option is correct
    is_B_correct = (chosen_option_config["B"] == required_B_stereochem)

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            error_messages.append(
                f"Constraint for starting material A is not satisfied. "
                f"Reaction A produces an (R)-lactone, and since stereochemistry is retained, "
                f"starting material A must be (R). The chosen answer suggests A is ({chosen_option_config['A']})."
            )
        if not is_B_correct:
            error_messages.append(
                f"Constraint for starting material B is not satisfied. "
                f"Reaction B produces an (S)-lactone, and since stereochemistry is retained, "
                f"starting material B must be (S). The chosen answer suggests B is ({chosen_option_config['B']})."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check
result = check_chemistry_answer()
print(result)