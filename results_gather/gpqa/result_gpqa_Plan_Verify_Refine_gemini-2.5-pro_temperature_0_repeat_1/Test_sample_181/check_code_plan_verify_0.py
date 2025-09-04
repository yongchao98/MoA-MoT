def check_answer_correctness():
    """
    Checks the correctness of the selected answer for the Mott-Gurney law question.
    """

    # Define the fundamental assumptions for the Mott-Gurney law's validity.
    # "no_injection_barrier" is a direct consequence and synonym for an ideal Ohmic contact in this context.
    mott_gurney_assumptions = {
        "carrier_model": "single-carrier",
        "trap_condition": "trap-free",
        "contact_condition": "no_injection_barrier", # or "ohmic"
        "dominant_current": "drift",
        "negligible_current": "diffusion"
    }

    # Parse the statements from each option into a structured format.
    options = {
        "A": {
            "carrier_model": "single-carrier",
            "trap_condition": "trap-free",
            "contact_condition": "ohmic",
            "negligible_current": "drift" # Contradicts the drift-dominated assumption.
        },
        "B": {
            "carrier_model": "single-carrier",
            "contact_condition": "schottky", # A Schottky contact implies an injection barrier.
            "negligible_current": "diffusion"
        },
        "C": {
            "carrier_model": "single-carrier",
            "trap_condition": "trap-free",
            "contact_condition": "no_injection_barrier",
            "negligible_current": "diffusion" # This implies drift is dominant.
        },
        "D": {
            "carrier_model": "two-carrier", # Contradicts the single-carrier assumption.
            "contact_condition": "ohmic",
            "negligible_current": "diffusion"
        }
    }

    llm_answer = "C"
    selected_option = options[llm_answer]

    # --- Verification Logic ---
    # 1. Check carrier model
    if selected_option.get("carrier_model") != mott_gurney_assumptions["carrier_model"]:
        return f"Incorrect: The answer states the device is '{selected_option.get('carrier_model')}', but the Mott-Gurney law is derived for a '{mott_gurney_assumptions['carrier_model']}' device."

    # 2. Check trap condition
    if selected_option.get("trap_condition") != mott_gurney_assumptions["trap_condition"]:
        return f"Incorrect: The answer does not specify a '{mott_gurney_assumptions['trap_condition']}' solid, which is a key assumption."

    # 3. Check contact condition
    contact_ok = (selected_option.get("contact_condition") == "ohmic" or
                  selected_option.get("contact_condition") == "no_injection_barrier")
    if not contact_ok:
        return f"Incorrect: The answer specifies a '{selected_option.get('contact_condition')}' contact. The Mott-Gurney law requires an Ohmic contact (no injection barrier) to ensure the current is space-charge limited, not injection-limited."

    # 4. Check current mechanism
    if selected_option.get("negligible_current") != mott_gurney_assumptions["negligible_current"]:
        return f"Incorrect: The answer states that '{selected_option.get('negligible_current')}' current is negligible. The Mott-Gurney law assumes '{mott_gurney_assumptions['dominant_current']}' current is dominant and '{mott_gurney_assumptions['negligible_current']}' current is negligible."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_answer_correctness()
print(result)