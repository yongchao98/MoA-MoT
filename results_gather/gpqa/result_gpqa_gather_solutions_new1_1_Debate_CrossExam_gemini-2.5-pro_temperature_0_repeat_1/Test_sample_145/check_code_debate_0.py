def check_diels_alder_product():
    """
    Checks the correctness of the answer for the Diels-Alder reaction of
    5-fluorocyclopenta-1,3-diene with maleic anhydride.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'D'

    # 1. Define the options and their structural features based on IUPAC nomenclature.
    options = {
        'A': {'stereochem': '(3aR,4R,7S,7aS,8r)'},
        'B': {'stereochem': '(3aR,4R,7S,7aS,8s)'},
        'C': {'stereochem': '(3aR,4S,7R,7aS,8s)'},
        'D': {'stereochem': '(3aR,4S,7R,7aS,8r)'}
    }

    # Interpret the stereochemical descriptors into structural properties.
    for opt, data in options.items():
        # Rule 1: Endo vs. Exo addition type.
        # (3aR,4S,7R,7aS) corresponds to the endo adduct.
        # (3aR,4R,7S,7aS) corresponds to the exo adduct.
        if '4S,7R' in data['stereochem']:
            data['addition_type'] = 'endo'
        elif '4R,7S' in data['stereochem']:
            data['addition_type'] = 'exo'
        else:
            data['addition_type'] = 'unknown'

        # Rule 2: Syn vs. Anti product type.
        # '8r' corresponds to the anti-product (F is opposite the anhydride ring).
        # '8s' corresponds to the syn-product (F is on the same side as the anhydride ring).
        if '8r' in data['stereochem']:
            data['product_type'] = 'anti'
        elif '8s' in data['stereochem']:
            data['product_type'] = 'syn'
        else:
            data['product_type'] = 'unknown'

    # 2. Apply chemical principles to determine the major product.

    # Principle 1: The Alder Endo Rule.
    # The major product is the kinetically favored 'endo' adduct.
    # We filter out the 'exo' products.
    endo_candidates = {opt: data for opt, data in options.items() if data['addition_type'] == 'endo'}

    if llm_final_answer not in endo_candidates:
        return (f"Incorrect. The provided answer '{llm_final_answer}' corresponds to an 'exo' adduct. "
                f"The major product of a Diels-Alder reaction is typically the 'endo' adduct due to kinetic control (the endo rule).")

    # Principle 2: Facial Selectivity.
    # For a 5-fluoro substituent, electronic effects favor 'syn'-facial attack,
    # which results in the 'anti'-product.
    # We filter the remaining candidates for the 'anti' product type.
    final_candidates = {opt: data for opt, data in endo_candidates.items() if data['product_type'] == 'anti'}

    if not final_candidates:
        return "Checking Error: No product satisfies both the endo rule and the facial selectivity rule."
    
    if len(final_candidates) > 1:
        return "Checking Error: Multiple products satisfy all chemical principles."

    correct_answer = list(final_candidates.keys())[0]

    # 3. Compare the derived correct answer with the LLM's answer.
    if llm_final_answer == correct_answer:
        return "Correct"
    else:
        # If the answer is wrong, provide a specific reason.
        if options[llm_final_answer]['addition_type'] != 'endo':
             # This case is already handled above, but included for completeness.
             return f"Incorrect. The answer {llm_final_answer} is an 'exo' adduct, but the major product is the 'endo' adduct."
        if options[llm_final_answer]['product_type'] != 'anti':
             return (f"Incorrect. The answer {llm_final_answer} is an 'endo, syn' product. "
                     f"While 'endo' is correct, the major product is the 'endo, anti' product. This results from the "
                     f"electronically favored 'syn'-facial attack, which leads to an 'anti' orientation of the fluorine "
                     f"substituent relative to the anhydride ring.")
        
        # Generic fallback
        return f"Incorrect. The predicted correct answer is {correct_answer}, but the provided answer is {llm_final_answer}."

# Run the check
result = check_diels_alder_product()
print(result)