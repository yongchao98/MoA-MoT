def check_chemistry_answer():
    """
    Checks the correctness of the answer to the stereochemistry problem.
    """
    # --- Problem Constraints ---
    # Reaction 1: A + LiBH4 -> (R)-lactone
    # Reaction 2: B + BH3   -> (S)-lactone
    required_product_config_A = 'R'
    required_product_config_B = 'S'

    # --- LLM's Answer ---
    llm_answer = "C"

    # --- Available Options ---
    # Each option defines the stereochemistry for starting materials A and B.
    options = {
        "A": {"A_config": "S", "B_config": "S"},
        "B": {"A_config": "R", "B_config": "R"},
        "C": {"A_config": "R", "B_config": "S"},
        "D": {"A_config": "S", "B_config": "R"}
    }

    if llm_answer not in options:
        return f"Invalid answer format. The answer '{llm_answer}' is not one of the valid options (A, B, C, D)."

    # Get the proposed starting materials from the LLM's answer
    proposed_configs = options[llm_answer]
    proposed_A_config = proposed_configs["A_config"]
    proposed_B_config = proposed_configs["B_config"]

    # --- Verification Step 1: Check Reaction A ---
    # Chemical Rule: The reaction with LiBH4 retains the stereochemistry of the starting material.
    # To get an (R) product, the starting material A must be (R).
    if proposed_A_config != required_product_config_A:
        return (f"Incorrect. The answer '{llm_answer}' is wrong for Reaction 1.\n"
                f"Reason: The question states that A + LiBH4 yields an (R)-lactone. "
                f"Since this reaction retains stereochemistry, the starting material A must have an (R) configuration. "
                f"The proposed answer suggests A is ({proposed_A_config}), which is incorrect.")

    # --- Verification Step 2: Check Reaction B ---
    # Chemical Rule: The reaction with BH3 also retains the stereochemistry of the starting material.
    # To get an (S) product, the starting material B must be (S).
    if proposed_B_config != required_product_config_B:
        return (f"Incorrect. The answer '{llm_answer}' is wrong for Reaction 2.\n"
                f"Reason: The question states that B + BH3 yields an (S)-lactone. "
                f"Since this reaction retains stereochemistry, the starting material B must have an (S) configuration. "
                f"The proposed answer suggests B is ({proposed_B_config}), which is incorrect.")

    # --- Conclusion ---
    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)