def check_diels_alder_product_stereochemistry():
    """
    Checks the correctness of the provided answer for the major product of the
    Diels-Alder reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The function codifies the following chemical principles:
    1. The Alder-Endo rule for kinetic control.
    2. The IUPAC stereodescriptors for endo/exo adducts.
    3. The syn-directive effect of electronegative substituents.
    4. The geometric outcome of a syn-facial attack.
    5. The Cahn-Ingold-Prelog (CIP) assignment for the C8 stereocenter.
    """
    
    # The final answer provided by the LLM being checked.
    provided_answer = "B"

    # Define the properties of each option based on its IUPAC name.
    # This requires correctly identifying the endo/exo skeleton and the relative
    # position (syn/anti) of the fluorine substituent.
    # - Endo skeleton: (3aR,4S,7R,7aS)
    # - Exo skeleton: (3aR,4R,7S,7aS)
    # - 8r descriptor: 'anti' position (F is opposite the anhydride ring)
    # - 8s descriptor: 'syn' position (F is on the same side as the anhydride ring)
    options_properties = {
        'A': {'skeleton': 'endo', 'substituent_pos': 'syn'},
        'B': {'skeleton': 'endo', 'substituent_pos': 'anti'},
        'C': {'skeleton': 'exo', 'substituent_pos': 'anti'},
        'D': {'skeleton': 'exo', 'substituent_pos': 'syn'}
    }

    # --- Step 1: Analyze Endo/Exo Selectivity ---
    # The Alder-Endo rule states the major kinetic product is the 'endo' adduct.
    major_skeleton = 'endo'
    
    possible_options = {opt: props for opt, props in options_properties.items() 
                        if props['skeleton'] == major_skeleton}

    # Check if the provided answer satisfies this first principle.
    if provided_answer not in possible_options:
        return (f"Incorrect. The provided answer '{provided_answer}' corresponds to an "
                f"'{options_properties[provided_answer]['skeleton']}' adduct. The major product must be an "
                f"'{major_skeleton}' adduct according to the Alder-Endo rule, which limits the "
                f"possibilities to {list(possible_options.keys())}.")

    # --- Step 2: Analyze Facial Selectivity and Resulting Geometry ---
    # For a 5-fluoro substituent, electronic effects favor 'syn'-facial attack.
    # A 'syn'-facial attack results in an 'anti'-adduct, where the fluorine on the C8 bridge
    # is on the opposite side of the bicyclic system from the anhydride ring.
    major_substituent_position = 'anti'

    # --- Step 3: Determine the correct option based on all principles ---
    final_correct_option = None
    for opt, props in possible_options.items():
        if props['substituent_pos'] == major_substituent_position:
            final_correct_option = opt
            break
            
    # --- Step 4: Compare the derived correct answer with the provided answer ---
    if provided_answer == final_correct_option:
        return "Correct"
    else:
        reason = (
            f"Incorrect. The provided answer is '{provided_answer}', but a rigorous analysis points to '{final_correct_option}'.\n"
            "Here is the correct line of reasoning:\n"
            f"1. **Endo/Exo Selectivity**: The Alder-Endo rule dictates that the major product is the '{major_skeleton}' adduct. This eliminates options C and D.\n"
            f"2. **Facial Selectivity**: For a 5-fluoro substituent, electronic effects favor 'syn'-facial attack.\n"
            f"3. **Product Geometry**: A 'syn'-facial attack leads to an '{major_substituent_position}' adduct, where the fluorine is on the opposite side of the anhydride ring.\n"
            f"4. **Conclusion**: The major product is the '{major_skeleton}, {major_substituent_position}' adduct. This corresponds to option '{final_correct_option}'. "
            f"The provided answer '{provided_answer}' corresponds to the '{options_properties[provided_answer]['skeleton']}, {options_properties[provided_answer]['substituent_pos']}' adduct, which is a minor product."
        )
        return reason

# Execute the check
result = check_diels_alder_product_stereochemistry()
print(result)