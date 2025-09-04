import re

def check_correctness():
    """
    This function checks the correctness of the provided answer about the Mott-Gurney equation.
    It defines the known physical constraints for the equation's validity and checks which option satisfies all of them.
    """

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    # Define the options from the question as a dictionary of their properties.
    # This allows for systematic checking of each condition.
    options = {
        'A': {
            "text": "The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
            "carrier_type": "single-carrier",
            "contact": "Schottky contact",
            "negligible_current": "diffusion current",
            "material": "not specified"
        },
        'B': {
            "text": "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.",
            "carrier_type": "single-carrier",
            "contact": "no carrier injection barrier",
            "negligible_current": "diffusion current",
            "material": "trap-free"
        },
        'C': {
            "text": "The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
            "carrier_type": "single-carrier",
            "contact": "Ohmic contact",
            "negligible_current": "drift current",
            "material": "trap-free"
        },
        'D': {
            "text": "The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
            "carrier_type": "two-carrier",
            "contact": "Ohmic contact",
            "negligible_current": "diffusion current",
            "material": "not specified"
        }
    }

    # Define the ground truth conditions for the Mott-Gurney equation's validity.
    # These are based on fundamental semiconductor physics.
    correct_conditions = {
        "carrier_type": "single-carrier",
        "contact": ["Ohmic contact", "no carrier injection barrier"], # Both phrases are valid for an ideal injecting contact.
        "negligible_current": "diffusion current", # The current is drift-dominated.
        "material": "trap-free"
    }

    # Find the option that satisfies all correct conditions.
    truly_correct_option_key = None
    for key, properties in options.items():
        # Check each condition
        cond1 = properties["carrier_type"] == correct_conditions["carrier_type"]
        cond2 = properties["contact"] in correct_conditions["contact"]
        cond3 = properties["negligible_current"] == correct_conditions["negligible_current"]
        cond4 = properties["material"] == correct_conditions["material"]

        if all([cond1, cond2, cond3, cond4]):
            truly_correct_option_key = key
            break

    # If no option is found to be correct, there's an issue with the question or the logic.
    if truly_correct_option_key is None:
        return "Error: The checking code could not find any option that satisfies all the physical constraints."

    # Compare the LLM's answer with the one derived from the physical constraints.
    if llm_answer == truly_correct_option_key:
        return "Correct"
    else:
        # If the LLM's answer is wrong, explain why.
        # First, explain why the LLM's choice is incorrect.
        llm_choice_properties = options[llm_answer]
        reasons = []
        if llm_choice_properties["carrier_type"] != correct_conditions["carrier_type"]:
            reasons.append(f"it specifies a '{llm_choice_properties['carrier_type']}' device, but the model is for a '{correct_conditions['carrier_type']}' device.")
        if llm_choice_properties["contact"] not in correct_conditions["contact"]:
            reasons.append(f"it specifies a '{llm_choice_properties['contact']}', but an ideal injecting contact ('Ohmic' or 'no injection barrier') is required.")
        if llm_choice_properties["negligible_current"] != correct_conditions["negligible_current"]:
            reasons.append(f"it states that '{llm_choice_properties['negligible_current']}' is negligible, but it should be '{correct_conditions['negligible_current']}'. The SCLC is a drift-dominated current.")
        if llm_choice_properties["material"] != correct_conditions["material"]:
            reasons.append(f"it does not specify a '{correct_conditions['material']}' material, which is a key assumption.")
        
        reason_string = " and ".join(reasons)
        return f"Incorrect. The provided answer is '{llm_answer}', but the correct answer is '{truly_correct_option_key}'. The answer '{llm_answer}' is wrong because {reason_string}"

# Execute the check and print the result.
result = check_correctness()
print(result)