def check_answer_correctness():
    """
    Checks if the provided answer correctly identifies the conditions for the
    validity of the Mott-Gurney equation.
    """
    # The final answer from the LLM analysis to be verified.
    llm_answer = "D"

    # The options as presented in the final analysis section of the prompt.
    options = {
        "A": "The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
        "B": "The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
        "C": "The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
        "D": "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current."
    }

    # Define the conditions for the correct answer.
    # A correct statement must not contain any 'forbidden' terms and must contain all 'required' terms to be complete.
    
    # Terms that are explicitly incorrect.
    forbidden_terms = {
        "two-carrier": "The law is for a single-carrier device.",
        "Schottky contact": "A Schottky contact has an injection barrier, but the law requires an Ohmic contact (no barrier).",
        "negligible drift current": "The current is drift-dominated; it is the diffusion current that is negligible."
    }

    # Terms that are required for a complete and correct statement.
    required_terms = [
        "trap-free",
        "single-carrier",
        "negligible diffusion"
    ]
    
    # The correct statement must also specify an Ohmic contact.
    required_contact_terms = [
        "Ohmic contact",
        "no carrier injection barrier"
    ]

    analysis_log = {}
    identified_correct_option = None

    for key, text in options.items():
        violations = []
        
        # 1. Check for any forbidden (contradictory) terms.
        for term, reason in forbidden_terms.items():
            if term in text:
                violations.append(reason)
        
        if violations:
            analysis_log[key] = (False, violations)
            continue

        # 2. If no contradictions, check if it's a complete answer.
        is_complete = True
        missing_conditions = []

        for term in required_terms:
            if term not in text:
                is_complete = False
                missing_conditions.append(f"is missing the required condition '{term}'")
        
        if not any(term in text for term in required_contact_terms):
            is_complete = False
            missing_conditions.append("is missing the required contact condition (Ohmic or no barrier)")

        if is_complete:
            # This option has no violations and is complete. It's the correct one.
            identified_correct_option = key
            analysis_log[key] = (True, [])
        else:
            # This option is not strictly wrong, but it is incomplete.
            analysis_log[key] = (False, missing_conditions)

    # 3. Compare the identified correct option with the LLM's answer.
    if identified_correct_option is None:
        return "Error: The checker could not identify a single, complete, and correct option among the choices."

    if llm_answer == identified_correct_option:
        return "Correct"
    else:
        is_llm_choice_correct, reasons = analysis_log[llm_answer]
        reason_str = " and ".join(reasons)
        return (f"Incorrect. The provided answer is '{llm_answer}', but the correct answer is '{identified_correct_option}'. "
                f"Option '{llm_answer}' is incorrect because: {reason_str}.")

# Run the checker and print the result.
result = check_answer_correctness()
print(result)