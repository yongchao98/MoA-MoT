def check_correctness():
    """
    Checks the correctness of the answer by verifying biological constraints.

    The riddle describes the secretory pathway for a protein.
    1.  Start: Protein synthesis begins on a free ribosome in the cytosol.
    2.  Pathway: The protein is targeted to the "rough" ER, where it gets "sugar" (glycosylated).
    3.  End: A common final destination for this pathway is secretion into the extracellular space.
    """

    # The options represent the start and end of the protein's journey.
    options = {
        "A": {"start": "ribosome", "end": "proteasome"},
        "B": {"start": "Golgi", "end": "mitochondrion"},
        "C": {"start": "membrane", "end": "nucleus"},
        "D": {"start": "cytosol", "end": "extracellular space"}
    }

    # The answer provided by the LLM.
    llm_answer_key = "D"

    # --- Constraint 1: The journey must start where protein synthesis begins for this pathway. ---
    # The SRP binds the nascent chain on a ribosome in the cytosol.
    valid_start_locations = ["cytosol", "ribosome"]
    
    # --- Constraint 2: The destination must be a valid endpoint for the secretory pathway. ---
    # The clues ("rough" ER, "sugar") point to the secretory pathway.
    # The proteasome is for degradation. Mitochondria and the nucleus have their own import pathways.
    # The extracellular space is a primary destination for secreted proteins.
    valid_secretory_destinations = ["extracellular space", "plasma membrane", "lysosome"]

    # Check the provided answer
    chosen_option = options.get(llm_answer_key)
    
    if not chosen_option:
        return f"Invalid answer key '{llm_answer_key}'. The key is not one of the options A, B, C, D."

    # Check start location
    if chosen_option["start"] not in valid_start_locations:
        return (f"Incorrect. The answer {llm_answer_key} ('{chosen_option['start']}' to '{chosen_option['end']}') is wrong. "
                f"Constraint Violated: The process described begins in the cytosol (on a ribosome), but the answer's start location is '{chosen_option['start']}'.")

    # Check end location
    if chosen_option["end"] not in valid_secretory_destinations:
        return (f"Incorrect. The answer {llm_answer_key} ('{chosen_option['start']}' to '{chosen_option['end']}') is wrong. "
                f"Constraint Violated: The clues ('rough' ER, 'sugar') indicate the secretory pathway. "
                f"'{chosen_option['end']}' is not a valid destination for this pathway. The proteasome (A), mitochondrion (B), and nucleus (C) all use different targeting systems.")

    # Verify that no other option is also correct
    for key, option_details in options.items():
        if key == llm_answer_key:
            continue
        if option_details["start"] in valid_start_locations and option_details["end"] in valid_secretory_destinations:
            return (f"Incorrect. While the provided answer {llm_answer_key} is plausible, option {key} also satisfies the constraints. "
                    f"The question may be ambiguous or the answer is not uniquely correct.")

    # If all checks pass and the answer is unique
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)