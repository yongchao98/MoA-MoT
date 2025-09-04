def check_diels_alder_stereochemistry():
    """
    This function checks the correctness of the answer to a Diels-Alder stereochemistry question.

    The logic follows established principles of pericyclic reactions:
    1.  The reaction is a Diels-Alder, which is kinetically controlled.
    2.  The 'endo rule' applies, favoring the endo adduct due to secondary orbital overlap.
    3.  For substituted cyclic dienes, the dienophile attack is sterically directed. It approaches
        from the face 'anti' (opposite) to the substituent on the sp3 carbon.
    4.  An 'endo' approach combined with an 'anti' attack results in a 'syn' orientation of the
        substituent relative to the newly formed anhydride bridge in the product.
    5.  The final step is to map this 'endo-syn' product to the correct IUPAC name among the options.
    """

    # The user's provided answer from the LLM
    llm_answer_option = "B"
    llm_reasoning = "The major product ... is the endo-syn adduct. ... results from the kinetically favored endo approach, where the dienophile attacks from the face anti to the fluorine substituent ... This 'anti-attack' geometry leads to a final product where the fluorine atom is syn to the anhydride bridge. ... The IUPAC name corresponding to the endo-syn isomer is (3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione, which is Option B."

    # --- Step 1: Define the expected outcome based on chemical principles ---

    # Principle 1: Endo/Exo selectivity
    # The kinetically favored product in a Diels-Alder reaction is the 'endo' adduct.
    favored_approach = "endo"

    # Principle 2: Syn/Anti attack selectivity
    # Steric hindrance from the C5-Fluorine favors attack 'anti' to the substituent.
    favored_attack = "anti"

    # Resulting product stereochemistry: endo approach + anti attack -> syn product
    major_product_type = "endo-syn"

    # --- Step 2: Map the options to their stereochemical descriptions ---
    # This mapping is based on IUPAC nomenclature for bicyclo[2.2.1]heptane systems.
    # - (...,4S,7R,...) denotes an 'endo' adduct.
    # - (...,4R,7S,...) denotes an 'exo' adduct.
    # - For the C8 substituent, '8r' corresponds to 'syn' and '8s' to 'anti' in this specific endo system.
    options_map = {
        "A": "endo-anti",
        "B": "endo-syn",
        "C": "exo-syn",  # From exo approach, anti attack
        "D": "exo-anti"  # From exo approach, syn attack
    }

    # --- Step 3: Determine the correct option letter ---
    correct_option = None
    for option, product_type in options_map.items():
        if product_type == major_product_type:
            correct_option = option
            break

    # --- Step 4: Validate the LLM's answer and reasoning ---
    if llm_answer_option != correct_option:
        return (f"Incorrect. The LLM's answer is {llm_answer_option}, but the correct option is {correct_option}. "
                f"The major product should be the '{major_product_type}' isomer, which corresponds to option {correct_option}.")

    # Check if the reasoning provided by the LLM is consistent with the principles
    if "endo-syn" not in llm_reasoning:
        return f"Incorrect. The LLM's answer {llm_answer_option} is correct, but the reasoning is flawed. It fails to identify the major product as the 'endo-syn' isomer."
    if "endo approach" not in llm_reasoning or "anti to the fluorine" not in llm_reasoning:
        return f"Incorrect. The LLM's answer {llm_answer_option} is correct, but the reasoning is incomplete. It does not correctly state that the 'endo-syn' product results from an 'endo' approach and an 'anti' attack."
    if options_map[llm_answer_option] != "endo-syn":
        return f"Incorrect. The LLM's reasoning correctly identifies the 'endo-syn' product as the major one, but incorrectly maps it to option {llm_answer_option}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_diels_alder_stereochemistry()
print(result)