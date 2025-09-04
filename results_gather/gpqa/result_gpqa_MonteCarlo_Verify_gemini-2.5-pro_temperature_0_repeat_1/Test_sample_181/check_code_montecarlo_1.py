def check_correctness():
    """
    This function checks the correctness of the given answer about the Mott-Gurney equation.
    It verifies the answer against the known physical assumptions required for the equation's validity.
    """
    # The given answer from the LLM
    llm_answer_key = 'D'

    # The options from the question
    options = {
        'A': "The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
        'B': "The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
        'C': "The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
        'D': "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current."
    }

    # Get the full text of the answer to be checked
    answer_text = options.get(llm_answer_key)

    if not answer_text:
        return f"Error: The answer key '{llm_answer_key}' is not a valid option."

    # Define the core assumptions for the Mott-Gurney Law.
    # The answer must satisfy all of these.
    # The contact condition has two valid phrasings.
    required_conditions = {
        "single-carrier": "The device must be single-carrier dominant.",
        "trap-free": "The material must be assumed to be trap-free.",
        "negligible diffusion current": "Diffusion current must be negligible compared to drift current.",
        "ohmic_contact_or_no_barrier": "There must be an Ohmic contact, which means no (or a negligible) carrier injection barrier."
    }

    # Define conditions that are explicitly incorrect.
    # The answer must not contain any of these.
    invalidating_conditions = {
        "two-carrier": "The Mott-Gurney law is for single-carrier, not two-carrier, devices.",
        "schottky contact": "A Schottky contact is rectifying and creates an injection barrier, violating the Ohmic contact assumption.",
        "negligible drift current": "The SCLC is a drift-dominated current. It is the diffusion current that is assumed to be negligible, not the drift current."
    }

    # Normalize the answer text for case-insensitive matching
    text_lower = answer_text.lower()
    
    # 1. Check for any explicitly incorrect statements
    for invalid_term, reason in invalidating_conditions.items():
        if invalid_term in text_lower:
            return f"Incorrect. The answer contains the invalid condition '{invalid_term}'. Reason: {reason}"

    # 2. Check for all required conditions
    # Check for 'single-carrier'
    if 'single-carrier' not in text_lower:
        return f"Incorrect. The answer is missing the required condition: 'single-carrier'. Reason: {required_conditions['single-carrier']}"
        
    # Check for 'trap-free'
    if 'trap-free' not in text_lower:
        return f"Incorrect. The answer is missing the required condition: 'trap-free'. Reason: {required_conditions['trap-free']}"

    # Check for 'negligible diffusion current'
    if 'negligible diffusion current' not in text_lower:
        return f"Incorrect. The answer is missing the required condition: 'negligible diffusion current'. Reason: {required_conditions['negligible diffusion current']}"

    # Check for the contact condition (Ohmic or no barrier)
    if 'ohmic contact' not in text_lower and 'no carrier injection barrier' not in text_lower:
        return f"Incorrect. The answer does not correctly describe the contact type. Reason: {required_conditions['ohmic_contact_or_no_barrier']}"

    # If all checks pass, the answer is correct.
    # A quick mental check of other options confirms they are wrong:
    # A: Fails on "negligible drift current".
    # B: Fails on "two-carrier".
    # C: Fails on "schottky contact" and missing "trap-free".
    # D: Passes all checks.
    
    return "Correct"

# Run the check and print the result
print(check_correctness())