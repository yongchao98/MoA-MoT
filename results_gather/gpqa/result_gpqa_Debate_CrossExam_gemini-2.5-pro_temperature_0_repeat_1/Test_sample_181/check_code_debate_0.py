def check_answer_correctness():
    """
    Checks the correctness of the LLM's answer about the Mott-Gurney equation.

    The function encodes the physical principles of the Mott-Gurney law and
    evaluates each option against these principles to determine the correct choice.
    """
    llm_answer = "C"

    # 1. Define the ground truth: the fundamental assumptions for the ideal Mott-Gurney law.
    ground_truth_conditions = {
        "carrier_type": "single-carrier",
        "trap_status": "trap-free",
        "contact_type": "ohmic",  # An ideal Ohmic contact provides no barrier to carrier injection.
        "dominant_current": "drift", # This implies that the diffusion current is negligible.
    }

    # 2. Parse the statements from each multiple-choice option.
    # Note: "no carrier injection barrier" is functionally equivalent to an ideal "ohmic" contact for SCLC.
    options = {
        "A": {
            "statement": "The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
            "conditions": {
                "carrier_type": "single-carrier",
                "contact_type": "schottky",
                "negligible_current": "diffusion"
            }
        },
        "B": {
            "statement": "The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
            "conditions": {
                "carrier_type": "single-carrier",
                "trap_status": "trap-free",
                "contact_type": "ohmic",
                "negligible_current": "drift"
            }
        },
        "C": {
            "statement": "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.",
            "conditions": {
                "carrier_type": "single-carrier",
                "trap_status": "trap-free",
                "contact_type": "no_injection_barrier",
                "negligible_current": "diffusion"
            }
        },
        "D": {
            "statement": "The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
            "conditions": {
                "carrier_type": "two-carrier",
                "contact_type": "ohmic",
                "negligible_current": "diffusion"
            }
        }
    }

    # 3. Systematically check each option.
    determined_correct_option = None
    error_log = {}

    for option_key, option_data in options.items():
        is_correct = True
        reasons_for_failure = []
        conds = option_data["conditions"]

        # Check carrier type
        if conds.get("carrier_type") != ground_truth_conditions["carrier_type"]:
            is_correct = False
            reasons_for_failure.append(f"it assumes a '{conds.get('carrier_type')}' device, but the law is for a '{ground_truth_conditions['carrier_type']}' device.")

        # Check trap status (the ideal law is for trap-free)
        if "trap_status" in conds and conds.get("trap_status") != ground_truth_conditions["trap_status"]:
            is_correct = False
            reasons_for_failure.append(f"it incorrectly specifies trap status.")
        
        # Check contact type
        contact = conds.get("contact_type")
        if contact == "schottky":
            is_correct = False
            reasons_for_failure.append("it specifies a Schottky contact, which is a barrier to injection, whereas an Ohmic contact (no barrier) is required for SCLC.")
        elif contact not in ["ohmic", "no_injection_barrier"]:
            is_correct = False
            reasons_for_failure.append(f"it specifies an invalid contact type '{contact}'.")

        # Check dominant current mechanism
        if conds.get("negligible_current") == "drift":
            is_correct = False
            reasons_for_failure.append("it states drift current is negligible, which is incorrect. The SCLC is a drift-dominated current.")
        elif conds.get("negligible_current") != "diffusion":
            is_correct = False
            reasons_for_failure.append("it incorrectly identifies the negligible current component.")

        # The most complete and correct option should be chosen.
        # Option C correctly identifies all key aspects: single-carrier, trap-free, no injection barrier, and negligible diffusion.
        if is_correct:
            # This logic assumes only one option can be fully correct.
            determined_correct_option = option_key
        else:
            error_log[option_key] = " and ".join(reasons_for_failure)

    # 4. Compare the derived correct answer with the LLM's answer.
    if llm_answer == determined_correct_option:
        return "Correct"
    else:
        if llm_answer not in options:
            return f"The provided answer '{llm_answer}' is not a valid option."
        
        failure_reason = error_log.get(llm_answer, "it is incomplete or incorrect compared to the best option.")
        return (f"The provided answer '{llm_answer}' is incorrect. \n"
                f"Reason: Option {llm_answer} is wrong because {failure_reason}\n"
                f"The correct answer is '{determined_correct_option}'. It is the only option that correctly states all the necessary conditions: a trap-free, single-carrier device with an Ohmic contact (no injection barrier) where the current is drift-dominated (i.e., diffusion current is negligible).")

# Execute the check and print the result.
result = check_answer_correctness()
print(result)