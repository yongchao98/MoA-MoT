def check_chemistry_answer():
    """
    This function checks the correctness of the answer to the chemistry question.
    It codifies the chemical principles of selectivity and stereochemistry for the given reactions.
    """

    # --- Define the Problem Constraints ---
    # Reaction A: A + LiBH4 -> (R)-product
    # Reaction B: B + BH3 -> (S)-product
    problem_constraints = {
        "A": {"reagent": "LiBH4", "product_stereo": "R"},
        "B": {"reagent": "BH3", "product_stereo": "S"}
    }

    # --- Define Chemical Principles ---
    # The key principle is that the stereochemistry at the C3 chiral center is retained
    # because the reactions (reduction at C1/C5 and cyclization) do not affect it.
    # stereochemistry_change can be "retained", "inverted", or "racemized".
    stereochemistry_change = "retained"

    # --- Deduce Required Starting Materials based on Principles ---
    deduced_A_stereo = None
    if stereochemistry_change == "retained":
        deduced_A_stereo = problem_constraints["A"]["product_stereo"]
    elif stereochemistry_change == "inverted":
        deduced_A_stereo = "S" if problem_constraints["A"]["product_stereo"] == "R" else "R"

    deduced_B_stereo = None
    if stereochemistry_change == "retained":
        deduced_B_stereo = problem_constraints["B"]["product_stereo"]
    elif stereochemistry_change == "inverted":
        deduced_B_stereo = "S" if problem_constraints["B"]["product_stereo"] == "R" else "R"

    # --- Define the Options from the Question ---
    options = {
        "A": {"A_stereo": "S", "B_stereo": "R"},
        "B": {"A_stereo": "R", "B_stereo": "S"},
        "C": {"A_stereo": "S", "B_stereo": "S"},
        "D": {"A_stereo": "R", "B_stereo": "R"}
    }

    # --- The Answer Provided by the LLM ---
    llm_answer = "B"

    # --- Verification Logic ---
    if llm_answer not in options:
        return f"Invalid answer option '{llm_answer}'. The options are A, B, C, D."

    chosen_option = options[llm_answer]

    # Check starting material A
    if chosen_option["A_stereo"] != deduced_A_stereo:
        return (f"Incorrect for starting material A. The answer claims A is ({chosen_option['A_stereo']}). "
                f"However, Reaction A produces an (R)-product. Since the reaction retains stereochemistry, "
                f"the starting material A must be (R) to yield the (R)-product.")

    # Check starting material B
    if chosen_option["B_stereo"] != deduced_B_stereo:
        return (f"Incorrect for starting material B. The answer claims B is ({chosen_option['B_stereo']}). "
                f"However, Reaction B produces an (S)-product. Since the reaction retains stereochemistry, "
                f"the starting material B must be (S) to yield the (S)-product.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_chemistry_answer()
print(result)