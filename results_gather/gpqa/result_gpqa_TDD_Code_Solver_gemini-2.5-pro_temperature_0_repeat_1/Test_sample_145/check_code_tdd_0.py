def check_answer_correctness():
    """
    This function checks the correctness of the given LLM answer for the
    Diels-Alder reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The check is based on two fundamental stereochemical principles:
    1. The Endo Rule: The kinetic product of a Diels-Alder reaction involving a
       cyclic diene is the 'endo' adduct.
    2. Facial Selectivity: The dienophile will approach from the face 'anti' to
       the substituent on the diene's bridge to minimize steric hindrance.

    The major product must therefore be the 'endo-anti' adduct.
    """
    
    # The options provided in the question
    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # The answer provided by the LLM
    llm_answer_key = "D"

    # --- Define Stereochemical Rules based on IUPAC Nomenclature ---

    # Rule 1: Identify 'endo' products.
    # The 'endo' adduct has the anhydride ring syn to the C8 bridge. This corresponds
    # to the relative stereochemistry (...,4S,7R,...) in the given options.
    # The 'exo' adduct would be (...,4R,7S,...).
    endo_stereochem_pattern = "(4S,7R,"

    # Rule 2: Identify 'anti' products.
    # 'Anti' attack means the fluorine at C8 is on the opposite side of the molecule
    # from the anhydride ring. This is denoted by the '8r' stereodescriptor.
    # 'Syn' attack would be denoted by '8s'.
    anti_stereochem_pattern = ",8r)"

    # --- Verification ---

    # Check if the provided answer key is valid
    if llm_answer_key not in options:
        return f"Invalid Answer: The key '{llm_answer_key}' is not one of the options (A, B, C, D)."

    selected_option_name = options[llm_answer_key]

    # Check if the selected option satisfies the 'endo' rule
    is_endo = endo_stereochem_pattern in selected_option_name
    if not is_endo:
        return (f"Incorrect: The answer '{llm_answer_key}' is wrong. "
                f"The major product must be the 'endo' adduct, which is kinetically favored. "
                f"This corresponds to the stereochemistry '{endo_stereochem_pattern}', but the chosen option is 'exo'.")

    # Check if the selected option satisfies the 'anti' rule
    is_anti = anti_stereochem_pattern in selected_option_name
    if not is_anti:
        return (f"Incorrect: The answer '{llm_answer_key}' is wrong. "
                f"The major product must result from 'anti' attack to minimize steric hindrance. "
                f"This corresponds to the stereodescriptor '{anti_stereochem_pattern}', but the chosen option is 'syn' ('8s').")

    # Final check: Ensure no other option also meets the criteria
    major_product_candidates = []
    for key, name in options.items():
        if endo_stereochem_pattern in name and anti_stereochem_pattern in name:
            major_product_candidates.append(key)
    
    if len(major_product_candidates) > 1:
         return (f"Ambiguous Question: The rules for identifying the major product are satisfied by multiple options: "
                 f"{major_product_candidates}. The question may be ill-posed.")

    if llm_answer_key not in major_product_candidates:
        # This case should be caught by the earlier checks, but serves as a failsafe.
        return (f"Incorrect: The answer '{llm_answer_key}' is wrong. "
                f"The correct option that is both 'endo' and 'anti' is '{major_product_candidates[0]}'.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_answer_correctness()
print(result)