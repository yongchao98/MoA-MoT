import re

def check_mott_gurney_answer(final_answer_text: str):
    """
    Checks the correctness of the final answer regarding the Mott-Gurney equation's validity.

    The function codifies the known physical assumptions of the Mott-Gurney law and
    evaluates which of the provided options (A, B, C, D) correctly describes them.
    It then compares this result with the given final answer.

    Args:
        final_answer_text: The full text of the final answer, which must end with <<<X>>>.

    Returns:
        A string indicating "Correct" or a reason for the incorrectness.
    """
    # Step 1: Extract the final answer letter (e.g., 'D') from the response.
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<X>>> in the provided text."
    final_answer_letter = match.group(1)

    # Step 2: Define the ground truth conditions for the Mott-Gurney equation.
    # These are the established physical assumptions for the model.
    ground_truth = {
        "carrier_type": "single-carrier",
        "material_purity": "trap-free",
        "contact_type": "no injection barrier",  # Synonymous with an ideal Ohmic contact for SCLC
        "negligible_current": "diffusion",
        "dominant_current": "drift"
    }

    # Step 3: Define the properties of each option based on the question's text.
    # The text for these options is taken from the final analysis provided in the prompt.
    options = {
        'A': {
            "text": "The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.",
            "properties": {
                "carrier_type": "single-carrier",
                "material_purity": "trap-free",
                "contact_type": "no injection barrier", # Ohmic contact implies this
                "negligible_current": "drift"
            }
        },
        'B': {
            "text": "The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.",
            "properties": {
                "carrier_type": "single-carrier",
                "material_purity": "not specified",
                "contact_type": "Schottky contact", # This is a barrier
                "negligible_current": "diffusion"
            }
        },
        'C': {
            "text": "The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.",
            "properties": {
                "carrier_type": "two-carrier",
                "material_purity": "not specified",
                "contact_type": "no injection barrier", # Ohmic contact implies this
                "negligible_current": "diffusion"
            }
        },
        'D': {
            "text": "The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.",
            "properties": {
                "carrier_type": "single-carrier",
                "material_purity": "trap-free",
                "contact_type": "no injection barrier",
                "negligible_current": "diffusion"
            }
        }
    }

    # Step 4: Determine the truly correct option by comparing each option's properties to the ground truth.
    correct_option_letter = None
    for letter, data in options.items():
        props = data["properties"]
        
        # Check all conditions
        is_correct = (
            props.get("carrier_type") == ground_truth["carrier_type"] and
            props.get("material_purity") == ground_truth["material_purity"] and
            props.get("contact_type") == ground_truth["contact_type"] and
            props.get("negligible_current") == ground_truth["negligible_current"]
        )
        
        if is_correct:
            correct_option_letter = letter
            break # Assuming only one option is fully correct

    # Step 5: Compare the provided final answer with the determined correct option.
    if final_answer_letter == correct_option_letter:
        return "Correct"
    else:
        # Provide a specific reason why the chosen answer is wrong.
        if final_answer_letter not in options:
             return f"The provided answer '{final_answer_letter}' is not a valid option choice."

        selected_props = options[final_answer_letter]["properties"]
        reasons = []
        if selected_props.get("carrier_type") != ground_truth["carrier_type"]:
            reasons.append(f"it incorrectly specifies a '{selected_props.get('carrier_type')}' device instead of a '{ground_truth['carrier_type']}' device")
        if selected_props.get("contact_type") != ground_truth["contact_type"]:
            reasons.append(f"it incorrectly specifies a '{selected_props.get('contact_type')}' instead of a contact with '{ground_truth['contact_type']}'")
        if selected_props.get("negligible_current") != ground_truth["negligible_current"]:
            reasons.append(f"it incorrectly states that '{selected_props.get('negligible_current')}' current is negligible, but it should be '{ground_truth['negligible_current']}' current (meaning drift current is dominant)")
        if selected_props.get("material_purity") != ground_truth["material_purity"] and selected_props.get("material_purity") != "not specified":
             reasons.append(f"it fails to specify the material must be '{ground_truth['material_purity']}'")

        if not reasons:
            return f"The provided answer '{final_answer_letter}' is incorrect. The correct answer is '{correct_option_letter}' based on a full analysis of the conditions."

        return f"The provided answer '{final_answer_letter}' is incorrect because " + ", and ".join(reasons) + "."

# The final answer provided by the LLM to be checked
final_answer_from_llm = """
Here is a step-by-step analysis to determine the correct answer.

The question asks for the set of conditions under which the Mott-Gurney equation is valid. The equation describes the Space-Charge-Limited Current (SCLC) regime. To find the correct statement, we must analyze the fundamental assumptions made in the derivation of this equation.

**1. Analyze the Core Assumptions of the Mott-Gurney Law:**

*   **Assumption 1: Carrier Type.** The model assumes that current is transported by only one type of charge carrier (either electrons or holes). The presence of a second mobile carrier type would introduce recombination effects, which are not included in this simple model. Therefore, the device must be a **single-carrier device**.
*   **Assumption 2: Transport Mechanism.** The model describes current driven by the electric field. It assumes the applied voltage is high enough that the **drift current** is the dominant transport mechanism. The contribution from **diffusion current** (due to carrier concentration gradients) is considered negligible.
*   **Assumption 3: Injection Contact.** For a space charge to build up and limit the current, the contact must be able to supply carriers into the device at a very high rate, essentially providing an unlimited reservoir of charge. This type of contact is called an **Ohmic contact**, which is characterized by having **no carrier injection barrier**. A contact with a significant barrier (like a Schottky contact) would limit the current by the injection process itself, not by the space charge in the bulk.
*   **Assumption 4: Material Purity.** The simplest form of the Mott-Gurney law assumes an ideal material with no localized energy states (traps) within the bandgap. If traps were present, they would capture some of the injected carriers, altering the current-voltage relationship. Therefore, the material must be **trap-free**.

**2. Evaluate Each Option Against the Assumptions:**

*   **A) The Mott-Gurney equation is valid for a trap-free single-carrier device with an Ohmic contact and negligible drift current.**
    *   This is **incorrect**. The Mott-Gurney equation describes a current that is *dominated* by drift. Stating that the drift current is negligible is a direct contradiction of the model.

*   **B) The Mott-Gurney equation is valid for a single-carrier device with a Schottky contact and negligible diffusion current.**
    *   This is **incorrect**. A Schottky contact has an energy barrier that impedes carrier injection. The Mott-Gurney law requires an Ohmic contact (no injection barrier) to ensure the current is limited by space charge, not by injection.

*   **C) The Mott-Gurney equation is valid for a two-carrier device with an Ohmic contact and negligible diffusion current.**
    *   This is **incorrect**. The Mott-Gurney law is a **single-carrier** model.

*   **D) The Mott-Gurney equation is valid for a trap-free single-carrier device with no carrier injection barrier and negligible diffusion current.**
    *   This is **correct**. This statement accurately combines all the necessary conditions:
        *   `trap-free`: Matches Assumption 4.
        *   `single-carrier device`: Matches Assumption 1.
        *   `no carrier injection barrier`: Matches Assumption 3.
        *   `negligible diffusion current`: Matches Assumption 2.

**Conclusion:**

By systematically checking each option against the fundamental physical assumptions required for the Mott-Gurney law, we can eliminate options A, B, and C. Option D is the only statement that correctly and completely describes the conditions for the equation's validity.

<<<D>>>
"""

# Run the check
result = check_mott_gurney_answer(final_answer_from_llm)
print(result)