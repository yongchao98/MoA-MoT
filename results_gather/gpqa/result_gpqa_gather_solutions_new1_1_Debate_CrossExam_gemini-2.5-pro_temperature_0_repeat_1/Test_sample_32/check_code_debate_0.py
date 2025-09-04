import re

def check_cycloaddition_answer():
    """
    Checks the correctness of the selected answer for the Diels-Alder reaction.
    This is done by verifying the product's structural features based on chemical principles.
    """

    # --- Problem Definition & Options ---
    question_details = {
        "diene": "2,5-dimethylthiophene",
        "target_product_type": "EXO"
    }

    options = {
        'A': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'C': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'D': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # The final answer provided by the LLM to be checked
    llm_answer = 'C'

    # --- Verification Rules based on Chemical Principles ---
    # For this reaction, the stereochemical descriptors for the EXO and ENDO products are known.
    stereochem_to_isomer_map = {
        "(3aR,4S,7R,7aS)": "EXO",
        "(3aR,4R,7S,7aS)": "ENDO"
    }

    # --- Verification Process ---
    selected_option_name = options.get(llm_answer)
    if not selected_option_name:
        return f"Error: The provided answer '{llm_answer}' is not a valid option key (A, B, C, or D)."

    # 1. Check the bridge type constraint
    # The diene is thiophene, so the bridge must be sulfur ("epithio").
    if "epithio" not in selected_option_name:
        # Find the incorrect term used
        incorrect_term = "epoxy" if "epoxy" in selected_option_name else "an incorrect"
        return (f"Incorrect. The answer '{llm_answer}' is wrong because the diene is 2,5-dimethylthiophene, "
                f"which must form an 'epithio' (sulfur) bridge. The name provided describes {incorrect_term} bridge.")

    # 2. Check the stereochemistry constraint
    # Extract the stereodescriptor, e.g., "(3aR,4S,7R,7aS)"
    match = re.search(r'^\(.*\)', selected_option_name)
    if not match:
        return f"Incorrect. The name for answer '{llm_answer}' is missing the required stereochemical descriptor at the beginning."

    descriptor = match.group(0)
    isomer_type = stereochem_to_isomer_map.get(descriptor)

    if not isomer_type:
        return (f"Incorrect. The stereochemical descriptor '{descriptor}' in answer '{llm_answer}' does not correspond "
                f"to a known standard configuration for this product.")

    if isomer_type != question_details["target_product_type"]:
        return (f"Incorrect. The answer '{llm_answer}' describes the {isomer_type} product, "
                f"but the question specifically asks for the {question_details['target_product_type']} product.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_cycloaddition_answer()
print(result)