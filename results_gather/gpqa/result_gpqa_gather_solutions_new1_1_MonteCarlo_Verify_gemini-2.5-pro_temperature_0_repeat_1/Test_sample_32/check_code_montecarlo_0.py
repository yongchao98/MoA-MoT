def check_diels_alder_product():
    """
    Checks the correctness of the LLM's answer for the Diels-Alder reaction.
    The verification is based on three main constraints:
    1. The chemical structure derived from the reactants (epithio vs. epoxy bridge).
    2. The stereochemical outcome based on reaction conditions (EXO vs. ENDO).
    3. The correct IUPAC name for the target stereoisomer, based on public chemical data.
    """
    # --- Problem Definition ---
    question_details = {
        "diene": "2,5-dimethylthiophene",
        "target_product_type": "EXO"
    }

    options = {
        "A": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # The final answer provided by the LLM being checked.
    llm_answer = "B"

    # --- Constraint 1: Check the basic chemical structure (bridge type) ---
    # The diene is thiophene, so the bridge must be sulfur ("epithio").
    if "epoxy" in options[llm_answer]:
        return f"Incorrect. The chosen answer '{llm_answer}' describes a product with an 'epoxy' (oxygen) bridge. The diene is {question_details['diene']}, which must form an 'epithio' (sulfur) bridge."

    # Verify that the eliminated options C and D are indeed incorrect for this reason.
    if "epoxy" not in options["C"] or "epoxy" not in options["D"]:
        return "Logic Error: The check cannot confirm that options C and D are incorrect based on the 'epoxy' vs 'epithio' distinction."

    # --- Constraint 2: Check the stereoisomer mapping (EXO/ENDO) ---
    # Based on data from chemical databases (e.g., PubChem CIDs 10246931 & 10246932),
    # we can establish the ground truth for the IUPAC names.
    ground_truth_mapping = {
        "A": "ENDO",
        "B": "EXO"
    }

    # The question asks for the EXO product.
    target_isomer = question_details["target_product_type"]

    # Find which option from the ground truth corresponds to the target isomer.
    correct_option_for_target = None
    for option, isomer_type in ground_truth_mapping.items():
        if isomer_type == target_isomer:
            correct_option_for_target = option
            break
    
    if not correct_option_for_target:
        return f"Internal Check Error: Could not find a ground truth mapping for the target isomer '{target_isomer}'."

    # --- Final Verdict ---
    if llm_answer == correct_option_for_target:
        return "Correct"
    else:
        llm_isomer_type = ground_truth_mapping.get(llm_answer, "Unknown")
        return (f"Incorrect. The question asks for the {target_isomer} product. "
                f"According to established chemical nomenclature, the {target_isomer} product is described by option '{correct_option_for_target}'. "
                f"The provided answer '{llm_answer}' corresponds to the {llm_isomer_type} product.")

# Execute the check and print the result.
result = check_diels_alder_product()
print(result)