import re

def check_cycloaddition_answer():
    """
    Checks the correctness of the provided answer for the Diels-Alder reaction.
    The check is performed in three steps:
    1. Verify the basic chemical structure from the IUPAC name.
    2. Identify the target stereoisomer (EXO).
    3. Match the stereoisomer to the correct IUPAC stereochemical descriptor.
    """
    # --- Problem Definition ---
    question = {
        "diene": "2,5-dimethylthiophene",
        "dienophile": "Furan-2,5-dione",
        "condition": "Heat",
        "product_type": "EXO"
    }
    options = {
        "A": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }
    llm_answer = "B"

    # --- Step 1: Check the basic chemical structure from the name ---
    # The product must have a sulfur bridge ("epithio") and an isobenzofuran core.
    # "epoxy" or "benzo[c]thiophene" are incorrect.
    valid_options_by_name = []
    for option_key, name in options.items():
        is_correct_structure = "epithio" in name and "isobenzofuran" in name
        if is_correct_structure:
            valid_options_by_name.append(option_key)
        elif option_key == llm_answer:
            reason = f"Answer {llm_answer} is incorrect due to fundamental structural naming errors. "
            if "epoxy" in name:
                reason += "The name contains 'epoxy' (oxygen bridge), but the diene was thiophene, so the product must have an 'epithio' (sulfur) bridge."
            if "benzo[c]thiophene" in name:
                reason += "The base name is incorrect; it should be derived from the dienophile, resulting in 'isobenzofuran-1,3-dione'."
            return reason

    # After this check, only options B and D should remain.
    if set(valid_options_by_name) != {'B', 'D'}:
        return f"Error in checker logic: Name-based filtering resulted in {valid_options_by_name}, expected ['B', 'D']."

    # --- Step 2 & 3: Check stereochemistry ---
    # Heat favors the thermodynamic EXO product. We need to identify which name corresponds to EXO.
    # Based on established nomenclature rules for this system:
    # EXO adduct descriptor: (3aR,4S,7R,7aS)
    # ENDO adduct descriptor: (3aR,4R,7S,7aS)
    
    exo_descriptor_pattern = r"\(3aR,4S,7R,7aS\)"
    correct_option_key = None
    
    for option_key in valid_options_by_name:
        if re.search(exo_descriptor_pattern, options[option_key]):
            correct_option_key = option_key
            break
            
    if correct_option_key is None:
        return "Error in checker logic: Could not find the known EXO stereodescriptor among the options."

    # --- Final Verdict ---
    if llm_answer == correct_option_key:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option_key}. "
                f"The question asks for the EXO product, which is favored by heat. "
                f"The IUPAC name for the EXO adduct is {options[correct_option_key]}, which is option {correct_option_key}. "
                f"The selected answer, option {llm_answer}, describes the ENDO adduct.")

# Execute the check and print the result
result = check_cycloaddition_answer()
print(result)