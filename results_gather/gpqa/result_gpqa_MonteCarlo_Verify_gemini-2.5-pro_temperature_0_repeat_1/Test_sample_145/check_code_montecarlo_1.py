def check_diels_alder_answer():
    """
    Checks the correctness of the answer for the Diels-Alder reaction of
    5-fluorocyclopenta-1,3-diene and maleic anhydride.
    """
    # The provided answer from the LLM
    llm_answer = "B"

    # --- Step 1: Determine the expected major product based on chemical principles ---
    # Principle 1: The 'endo' rule favors the endo product.
    # Principle 2: Steric hindrance favors 'anti-attack' (dienophile approaches from the face opposite the fluorine).
    # Principle 3: 'Anti-attack' leads to a 'syn' product (fluorine is on the same side as the anhydride bridge).
    expected_major_product_type = "endo-syn"

    # --- Step 2: Map the IUPAC names in the options to their stereochemical descriptors ---
    # This mapping is based on Cahn-Ingold-Prelog rules for this specific bicyclic system.
    # endo core: (3aR,4S,7R,7aS)
    # exo core: (3aR,4R,7S,7aS)
    # syn fluorine (relative to anhydride bridge): 8r
    # anti fluorine (relative to anhydride bridge): 8s
    option_to_isomer_map = {
        "A": "endo-anti",
        "B": "endo-syn",
        "C": "exo-syn",
        "D": "exo-anti"
    }

    # --- Step 3: Verify the LLM's answer ---
    # Find the isomer type corresponding to the LLM's answer.
    llm_answer_isomer_type = option_to_isomer_map.get(llm_answer)

    if llm_answer_isomer_type is None:
        return f"Invalid option '{llm_answer}' provided."

    # Check if the LLM's answer matches the expected major product.
    if llm_answer_isomer_type == expected_major_product_type:
        return "Correct"
    else:
        # Find which option *should* have been correct.
        correct_option = None
        for option, isomer_type in option_to_isomer_map.items():
            if isomer_type == expected_major_product_type:
                correct_option = option
                break
        
        reason = (
            f"The answer is incorrect. The provided answer is {llm_answer}, which corresponds to the "
            f"'{llm_answer_isomer_type}' isomer. However, the reaction's major product is determined by "
            f"endo-selectivity and anti-attack, which results in the '{expected_major_product_type}' isomer. "
            f"The correct option corresponding to the '{expected_major_product_type}' isomer is {correct_option}."
        )
        return reason

# Run the check
result = check_diels_alder_answer()
print(result)