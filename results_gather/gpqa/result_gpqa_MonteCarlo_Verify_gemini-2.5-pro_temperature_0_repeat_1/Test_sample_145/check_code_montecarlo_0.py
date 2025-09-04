def check_diels_alder_product():
    """
    Checks the correctness of the LLM's answer by applying deterministic
    chemical principles of the Diels-Alder reaction.
    """

    # Principle 1: Determine the favored ring stereochemistry (Endo/Exo).
    # The Alder-endo rule predicts the 'endo' product is kinetically favored.
    favored_ring_form = "endo"

    # Principle 2: Determine the favored substituent orientation (Syn/Anti).
    # The dienophile attacks 'anti' to the C5-substituent (F) to minimize steric hindrance.
    # This 'anti-attack' leads to a 'syn' configuration in the final product, where the
    # F atom is on the opposite side of the molecule from the anhydride ring.
    favored_substituent_orientation = "syn"

    # The expected major product is the combination of these two features.
    expected_product_type = f"{favored_ring_form}-{favored_substituent_orientation}"

    # Provided options and the LLM's answer.
    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }
    llm_answer_choice = "A"

    # Map IUPAC names to stereochemical types based on Cahn-Ingold-Prelog rules.
    # This mapping is consistent with the one used in the provided answer.
    def get_product_type(name):
        if "(3aR,4S,7R,7aS" in name:
            ring = "endo"
        elif "(3aR,4R,7S,7aS" in name:
            ring = "exo"
        else:
            return "Unknown Ring"

        if ",8s)" in name:
            substituent = "anti"
        elif ",8r)" in name:
            substituent = "syn"
        else:
            return "Unknown Substituent"
        
        return f"{ring}-{substituent}"

    # Determine the stereochemical type of the LLM's chosen answer.
    llm_answer_type = get_product_type(options[llm_answer_choice])

    # Compare the expected product with the LLM's answer.
    if llm_answer_type == expected_product_type:
        return "Correct"
    else:
        correct_option = "Unknown"
        for option, name in options.items():
            if get_product_type(name) == expected_product_type:
                correct_option = option
                break
        
        reason = (
            f"The provided answer '{llm_answer_choice}' is incorrect.\n"
            f"1. **Correct Principles**: The major product is correctly identified as being the 'endo' adduct, resulting from a sterically favored 'anti-attack' of the dienophile relative to the fluorine atom.\n"
            f"2. **Incorrect Conclusion**: The answer incorrectly concludes that an 'anti-attack' leads to an 'anti' product. An 'anti-attack' (where the dienophile approaches from the side opposite the substituent) results in a 'syn' product configuration, where the substituent (F) and the new anhydride bridge are on opposite faces of the bicyclic system.\n"
            f"3. **Correct Product**: The major product should be the 'endo-syn' isomer.\n"
            f"4. **Correct Option**: The 'endo-syn' isomer corresponds to Option B, which has the stereodescriptors (3aR,4S,7R,7aS) for the endo core and '8r' for the syn fluorine.\n"
            f"5. **LLM's Error**: The LLM chose Option A, the 'endo-anti' isomer. This is a minor product."
        )
        return reason

# Execute the check and print the result.
result = check_diels_alder_product()
print(result)