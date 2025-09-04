def check_correctness():
    """
    Checks the correctness of the selected answer for the Mott-Gurney law question.

    The function codifies the known physical assumptions for the Mott-Gurney law
    and evaluates each option against them.
    """

    # Define the ground-truth conditions for the Mott-Gurney law's validity.
    required_conditions = {
        "carrier_type": "single-carrier",
        "material_condition": "trap-free",
        "contact_type": "no_injection_barrier",  # Synonymous with an ideal Ohmic contact for SCLC
        "current_mechanism": "negligible_diffusion" # Implies drift-dominated
    }

    # Define the properties of each option as presented in the final analysis.
    options = {
        "A": {
            "carrier_type": "single-carrier",
            "material_condition": "trap-free",
            "contact_type": "no_injection_barrier",
            "current_mechanism": "negligible_diffusion"
        },
        "B": {
            "carrier_type": "two-carrier",
            "material_condition": "not_specified",
            "contact_type": "ohmic",
            "current_mechanism": "negligible_diffusion"
        },
        "C": {
            "carrier_type": "single-carrier",
            "material_condition": "not_specified",
            "contact_type": "schottky",
            "current_mechanism": "negligible_diffusion"
        },
        "D": {
            "carrier_type": "single-carrier",
            "material_condition": "trap-free",
            "contact_type": "ohmic",
            "current_mechanism": "negligible_drift" # This is the key error in this option
        }
    }

    # The final answer provided by the LLM to be checked.
    llm_final_answer = "A"

    # Find the correct option based on the required conditions.
    correct_option_key = None
    for key, properties in options.items():
        # Check carrier type
        is_carrier_type_correct = (properties["carrier_type"] == required_conditions["carrier_type"])
        
        # Check material condition (must be trap-free)
        is_material_correct = (properties["material_condition"] == required_conditions["material_condition"])
        
        # Check contact type (must be Ohmic/no barrier, not Schottky)
        is_contact_correct = (properties["contact_type"] in ["ohmic", "no_injection_barrier"])
        
        # Check current mechanism (must be negligible diffusion, not negligible drift)
        is_current_mechanism_correct = (properties["current_mechanism"] == required_conditions["current_mechanism"])

        # An option is correct if and only if all conditions are met.
        if is_carrier_type_correct and is_material_correct and is_contact_correct and is_current_mechanism_correct:
            correct_option_key = key
            break # Assume only one correct answer

    # Verify the LLM's answer
    if llm_final_answer == correct_option_key:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness.
        if correct_option_key is None:
            return "Error in checking logic: No option satisfies all the required conditions."
        
        reasons = []
        props = options[llm_final_answer]
        if props["carrier_type"] != required_conditions["carrier_type"]:
            reasons.append(f"it describes a '{props['carrier_type']}' device instead of a '{required_conditions['carrier_type']}' device")
        if props["material_condition"] != required_conditions["material_condition"]:
            reasons.append(f"it fails to specify the material must be '{required_conditions['material_condition']}'")
        if props["contact_type"] not in ["ohmic", "no_injection_barrier"]:
            reasons.append(f"it specifies a '{props['contact_type']}' contact, which has an injection barrier")
        if props["current_mechanism"] != required_conditions["current_mechanism"]:
            reasons.append(f"it incorrectly states that '{props['current_mechanism'].split('_')[1]}' current is negligible, when it should be '{required_conditions['current_mechanism'].split('_')[1]}'")

        return (f"Incorrect. The provided answer '{llm_final_answer}' is wrong because {', and '.join(reasons)}. "
                f"The correct answer is '{correct_option_key}', which correctly states all conditions: single-carrier, trap-free, no injection barrier, and negligible diffusion current.")

# Run the check and print the result.
result = check_correctness()
print(result)