import re

def check_answer():
    """
    Checks the correctness of the LLM's answer by verifying chemical constraints.
    """
    # The options provided in the problem description
    options = {
        "A": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }

    # The final answer provided by the LLM being checked
    llm_answer_key = "B"
    llm_answer_text = options.get(llm_answer_key)

    if not llm_answer_text:
        return f"Invalid answer key '{llm_answer_key}'. Please provide a key from {list(options.keys())}."

    # --- Constraint 1: Check for correct connectivity and naming ---
    # The diene is thiophene, so the bridge must be sulfur ("epithio").
    # The dienophile is maleic anhydride, which forms an "isobenzofuran" core.
    required_name_parts = ["epithio", "isobenzofuran"]
    for part in required_name_parts:
        if part not in llm_answer_text:
            return (f"Incorrect. The answer '{llm_answer_key}' is wrong because its name is structurally incorrect. "
                    f"The product of 2,5-dimethylthiophene and maleic anhydride should contain the term '{part}', "
                    f"but the name is '{llm_answer_text}'.")

    # Check for incorrect parts, like "epoxy" (oxygen bridge)
    if "epoxy" in llm_answer_text:
        return (f"Incorrect. The answer '{llm_answer_key}' is wrong because it describes an 'epoxy' (oxygen) bridge. "
                f"The diene was thiophene, so the product must have an 'epithio' (sulfur) bridge.")

    # --- Constraint 2: Check for correct EXO stereochemistry ---
    # Based on rigorous Cahn-Ingold-Prelog analysis, the EXO adduct has a specific stereochemistry.
    # The reaction produces a racemic mixture, so we must check for either enantiomer.
    exo_stereochem_descriptors = [
        "(3aR,4S,7R,7aS)",  # One enantiomer
        "(3aS,4R,7S,7aR)"   # The other enantiomer
    ]
    
    # Extract the stereodescriptor from the name using regex
    match = re.match(r"(\([\w,]+\))", llm_answer_text)
    if not match:
        return f"Incorrect. The answer '{llm_answer_key}' does not contain a recognizable stereochemical descriptor."
        
    provided_descriptor = match.group(1)

    if provided_descriptor not in exo_stereochem_descriptors:
        # For completeness, let's identify the ENDO descriptor
        endo_stereochem_descriptors = ["(3aR,4R,7S,7aS)", "(3aS,4S,7R,7aR)"]
        if provided_descriptor in endo_stereochem_descriptors:
            return (f"Incorrect. The answer '{llm_answer_key}' is wrong because it describes the ENDO product. "
                    f"The question specifically asks for the EXO product, which has the stereochemistry "
                    f"'{exo_stereochem_descriptors[0]}' (or its enantiomer).")
        else:
            return (f"Incorrect. The answer '{llm_answer_key}' has an unrecognized stereochemical descriptor "
                    f"'{provided_descriptor}'.")

    # --- Final Verdict ---
    # If both connectivity and stereochemistry are correct
    return "Correct"

# Run the check
result = check_answer()
print(result)