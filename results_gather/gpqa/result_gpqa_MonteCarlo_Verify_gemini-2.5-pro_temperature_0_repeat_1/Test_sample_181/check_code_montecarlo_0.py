import collections

def check_mott_gurney_answer():
    """
    Checks the correctness of the provided answer for the Mott-Gurney equation conditions.

    The function defines the ground truth conditions for the Mott-Gurney law based on
    semiconductor physics and evaluates the given answer ('D') against these principles.
    If the answer satisfies all constraints, it returns "Correct". Otherwise, it
    returns the reasons why the answer is incorrect.
    """

    # The answer provided by the other LLM that we need to check.
    llm_answer = 'D'

    # Define the fundamental assumptions for the ideal Mott-Gurney Law.
    # This is the "ground truth" against which the answer will be checked.
    ground_truth_conditions = {
        "carrier_model": "single-carrier",
        "trap_status": "trap-free",
        "contact_type": "no injection barrier",  # This is the key property of an ideal Ohmic contact for SCLC.
        "dominant_transport": "drift",  # This implies that diffusion current is negligible.
    }

    # We parse the text of each option to extract its claimed conditions.
    # This allows for a systematic comparison with the ground truth.
    # We use a namedtuple for better readability.
    Option = collections.namedtuple('Option', ['text', 'conditions'])
    options_data = {
        'A': Option(
            text="The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
            conditions={
                "carrier_model": "single-carrier",
                "trap_status": "trap-free",
                "contact_type": "Ohmic contact",
                "dominant_transport": "diffusion", # FLAW: "negligible drift current" means diffusion is dominant.
            }
        ),
        'B': Option(
            text="The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
            conditions={
                "carrier_model": "two-carrier", # FLAW: The law is for single-carrier devices.
                "trap_status": "unspecified",
                "contact_type": "Ohmic contact",
                "dominant_transport": "drift",
            }
        ),
        'C': Option(
            text="The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
            conditions={
                "carrier_model": "single-carrier",
                "trap_status": "unspecified",
                "contact_type": "Schottky contact", # FLAW: A Schottky contact has an injection barrier.
                "dominant_transport": "drift",
            }
        ),
        'D': Option(
            text="The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.",
            conditions={
                "carrier_model": "single-carrier",
                "trap_status": "trap-free",
                "contact_type": "no injection barrier",
                "dominant_transport": "drift",
            }
        )
    }

    # --- Verification Logic ---
    if llm_answer not in options_data:
        return f"Error: The provided answer '{llm_answer}' is not a valid option (A, B, C, or D)."

    selected_option_conditions = options_data[llm_answer].conditions
    error_messages = []

    # 1. Check Carrier Model Constraint
    if selected_option_conditions.get("carrier_model") != ground_truth_conditions["carrier_model"]:
        error_messages.append(
            f"Constraint failure on carrier model: The law applies to a '{ground_truth_conditions['carrier_model']}' device, but the answer specifies '{selected_option_conditions.get('carrier_model')}'."
        )

    # 2. Check Trap Status Constraint
    if selected_option_conditions.get("trap_status") != ground_truth_conditions["trap_status"]:
        error_messages.append(
            f"Constraint failure on trap status: The ideal law assumes a '{ground_truth_conditions['trap_status']}' material, but the answer either does not specify this or specifies something else."
        )

    # 3. Check Contact Type Constraint
    # A Schottky contact is explicitly wrong. "Ohmic" is generally acceptable, but "no injection barrier" is more precise and correct.
    if selected_option_conditions.get("contact_type") == "Schottky contact":
        error_messages.append(
            "Constraint failure on contact type: The law requires a non-blocking contact (Ohmic/no injection barrier), but the answer specifies a 'Schottky contact' which has an injection barrier."
        )

    # 4. Check Dominant Transport Mechanism Constraint
    if selected_option_conditions.get("dominant_transport") != ground_truth_conditions["dominant_transport"]:
        error_messages.append(
            f"Constraint failure on transport mechanism: The law describes a '{ground_truth_conditions['dominant_transport']}-dominated' current, but the answer implies it is '{selected_option_conditions.get('dominant_transport')}-dominated'."
        )

    # --- Final Verdict ---
    if not error_messages:
        # The selected answer's conditions perfectly match all ground truth constraints.
        return "Correct"
    else:
        # The selected answer violates one or more physical constraints.
        reason = f"Incorrect. The answer '{llm_answer}' does not satisfy all constraints for the following reasons:\n"
        for i, msg in enumerate(error_messages, 1):
            reason += f"{i}. {msg}\n"
        return reason.strip()

# Execute the check and print the result.
result = check_mott_gurney_answer()
print(result)