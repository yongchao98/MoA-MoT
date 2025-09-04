def check_mott_gurney_validity():
    """
    This function checks the correctness of the provided answer 'D' for the question
    about the validity conditions of the Mott-Gurney equation. It does this by
    defining the known physical principles (ground truth) and verifying that the
    statement in option 'D' aligns with all of them.
    """
    
    # The final answer from the LLM to be verified.
    llm_answer = "D"

    # Ground Truth: The four key assumptions for the classic Mott-Gurney Law.
    correct_conditions = {
        "carrier_model": "single-carrier",
        "material_condition": "trap-free",
        "contact_type": "no carrier injection barrier",  # The defining feature of an ideal Ohmic contact for SCLC.
        "current_mechanism": "negligible diffusion"  # This implies the current is drift-dominated.
    }

    # The text of the option selected by the LLM.
    chosen_option_text = "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current."

    # A list to store reasons why the answer might be wrong.
    validation_errors = []

    # 1. Check for the correct carrier model.
    if correct_conditions["carrier_model"] not in chosen_option_text:
        validation_errors.append("The statement does not specify a 'single-carrier' device.")
    if "two-carrier" in chosen_option_text:
        validation_errors.append("The statement incorrectly refers to a 'two-carrier' device.")

    # 2. Check for the trap-free condition.
    if correct_conditions["material_condition"] not in chosen_option_text:
        validation_errors.append("The statement omits the essential 'trap-free' condition.")

    # 3. Check for the correct contact type.
    if correct_conditions["contact_type"] not in chosen_option_text:
        validation_errors.append("The statement omits the essential 'no carrier injection barrier' condition.")
    if "Schottky contact" in chosen_option_text:
        validation_errors.append("The statement incorrectly refers to a 'Schottky contact', which has an injection barrier.")

    # 4. Check for the correct current transport mechanism.
    if correct_conditions["current_mechanism"] not in chosen_option_text:
        validation_errors.append("The statement omits the 'negligible diffusion' current condition.")
    if "negligible drift" in chosen_option_text:
        validation_errors.append("The statement incorrectly claims 'negligible drift current', contradicting that SCLC is a drift-dominated current.")

    # Final verdict.
    if not validation_errors:
        return "Correct"
    else:
        error_summary = "; ".join(validation_errors)
        return f"Incorrect. The answer '{llm_answer}' is wrong because: {error_summary}"

# Execute the check and print the result.
result = check_mott_gurney_validity()
print(result)