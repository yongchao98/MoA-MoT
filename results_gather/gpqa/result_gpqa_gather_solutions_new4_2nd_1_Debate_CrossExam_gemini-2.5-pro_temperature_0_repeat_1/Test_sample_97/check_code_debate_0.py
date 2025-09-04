import re

def check_chemistry_answer():
    """
    Checks the correctness of the answer for a Ring-Opening Cross-Metathesis (ROCM) reaction.

    The function analyzes the reaction:
    A + methyleneruthenium compound + 1-propene ---> 1-(prop-1-en-1-yl)-2-vinylcyclopentane

    It evaluates the options based on established chemical principles for ROCM:
    1. The reaction is a Ring-Opening Cross-Metathesis (ROCM).
    2. The product has an intact cyclopentane core, so the starting material must contain one.
    3. The product is 1,2-disubstituted, meaning the ring that opened was fused to adjacent carbons of the cyclopentane core.
    4. The starting material must be a bicyclic alkene with a strained double bond suitable for ring-opening.
    5. The other options must be chemically inconsistent with the reaction outcome.
    """

    # --- Problem Definition ---
    question = "Identify the starting material, A, in the following reaction: A + a methyleneruthenium compound + 1-propene ---> 1-(prop-1-en-1-yl)-2-vinylcyclopentane"
    
    options = {
        "A": "2-methyl-3-methylenebicyclo[2.1.0]pentane",
        "B": "2-methylbicyclo[3.1.0]hex-2-ene",
        "C": "bicyclo[3.2.0]hept-6-ene",
        "D": "1,2-dimethylenecyclopentane"
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_text = "<<<C>>>"

    # --- Analysis of the Product and Reaction ---
    product_info = {
        "core_ring": "cyclopentane",
        "substitution_pattern": "1,2-disubstituted",
        "substituents": ["vinyl", "prop-1-en-1-yl"]
    }
    
    reaction_type = "Ring-Opening Cross-Metathesis (ROCM)"

    # --- Define Chemical Properties of Each Option ---
    # This simulates a chemist's analysis of each structure.
    option_properties = {
        "A": {
            "name": "2-methyl-3-methylenebicyclo[2.1.0]pentane",
            "is_bicyclic": True,
            "contains_cyclopentane_core": False, # Core is bicyclo[2.1.0]pentane
            "can_undergo_rocm_to_product": False
        },
        "B": {
            "name": "2-methylbicyclo[3.1.0]hex-2-ene",
            "is_bicyclic": True,
            "contains_cyclopentane_core": True,
            # The double bond is in the 5-membered ring. Opening it would destroy the core.
            "can_undergo_rocm_to_product": False 
        },
        "C": {
            "name": "bicyclo[3.2.0]hept-6-ene",
            "is_bicyclic": True,
            "contains_cyclopentane_core": True,
            # The double bond is in the strained 4-membered ring, perfect for ROCM.
            # Opening it creates a 1,2-disubstituted cyclopentane intermediate.
            "can_undergo_rocm_to_product": True
        },
        "D": {
            "name": "1,2-dimethylenecyclopentane",
            "is_bicyclic": False, # It's a diene, not bicyclic.
            "contains_cyclopentane_core": True,
            # Cannot undergo ring-opening.
            "can_undergo_rocm_to_product": False
        }
    }

    # --- Logic to Determine the Correct Answer ---
    correct_option = None
    reasons = []

    for key, props in option_properties.items():
        # Condition 1: Must be a bicyclic alkene to undergo ROCM.
        if not props["is_bicyclic"]:
            reasons.append(f"Option {key} ({props['name']}) is incorrect because it is not a bicyclic alkene and cannot undergo ring-opening.")
            continue
        
        # Condition 2: Must contain a cyclopentane core that remains intact.
        if not props["contains_cyclopentane_core"]:
            reasons.append(f"Option {key} ({props['name']}) is incorrect because it does not contain a cyclopentane ring system that can form the product's core.")
            continue

        # Condition 3: Must be able to form the product via ROCM. This means the double bond
        # must be in the strained ring, not the core ring.
        if props["can_undergo_rocm_to_product"]:
            if correct_option is None:
                correct_option = key
                reasons.append(f"Option {key} ({props['name']}) is the correct precursor. Its strained cyclobutene ring opens via ROCM to form a 1,2-disubstituted cyclopentane intermediate, which then reacts with 1-propene to yield the final product.")
            else:
                # This case handles if multiple options seem correct, which would indicate an issue.
                return "Error in analysis: Multiple options appear to be correct."
        else:
             reasons.append(f"Option {key} ({props['name']}) is incorrect. For B, ring-opening would destroy the cyclopentane core. For A, the core is wrong.")


    # --- Final Check ---
    # Extract the letter from the LLM's answer, e.g., 'C' from '<<<C>>>'
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format: {llm_answer_text}"
    
    llm_choice = match.group(1)

    if llm_choice == correct_option:
        return "Correct"
    else:
        failure_reason = f"The provided answer '{llm_choice}' is incorrect. The correct answer is '{correct_option}'.\n"
        # Find the specific reason for the correct choice.
        correct_reason = next((r for r in reasons if correct_option in r), "")
        failure_reason += correct_reason
        return failure_reason

# Execute the check and print the result
result = check_chemistry_answer()
print(result)