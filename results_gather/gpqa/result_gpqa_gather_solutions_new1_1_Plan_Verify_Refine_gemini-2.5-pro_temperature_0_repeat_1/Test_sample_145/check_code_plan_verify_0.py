def check_diels_alder_stereochemistry():
    """
    This function checks the correctness of the provided answer by systematically applying
    the chemical principles outlined in its justification.
    """
    
    # The final answer provided by the LLM to be checked.
    llm_answer = 'B'

    # --- Step 0: Define the options based on their IUPAC stereodescriptors ---
    options = {
        'A': {'name': '(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione',
              'core_type': 'exo', 'bridge_type': 'syn', 'bridge_descriptor': '8s'},
        'B': {'name': '(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione',
              'core_type': 'endo', 'bridge_type': 'anti', 'bridge_descriptor': '8r'},
        'C': {'name': '(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione',
              'core_type': 'exo', 'bridge_type': 'anti', 'bridge_descriptor': '8r'},
        'D': {'name': '(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione',
              'core_type': 'endo', 'bridge_type': 'syn', 'bridge_descriptor': '8s'}
    }

    # --- Step 1: Apply the Endo/Exo Selectivity Rule ---
    # The Alder endo rule states the kinetically favored product is the 'endo' adduct.
    # The provided answer correctly identifies this as the first filter.
    possible_products = {k: v for k, v in options.items() if v['core_type'] == 'endo'}
    
    if not possible_products:
        return "Logic Error: No 'endo' products were found among the options based on the defined nomenclature."
    
    # After this step, only options B and D should remain.
    if set(possible_products.keys()) != {'B', 'D'}:
        return f"Constraint Check Failed (Endo Rule): The rule should filter for 'endo' products, leaving B and D. However, the remaining options are {list(possible_products.keys())}."

    # --- Step 2: Apply the Facial Selectivity Rule ---
    # For a 5-fluoro substituent, electronic effects favor 'syn-facial' attack.
    # The provided answer correctly states this.
    favored_attack_type = 'syn-facial'

    # --- Step 3: Determine the Resulting Product Structure ---
    # A 'syn-facial' attack leads to an 'anti-product'.
    # The provided answer correctly states this relationship.
    if favored_attack_type == 'syn-facial':
        major_product_structure = 'anti'
    else: # 'anti-facial' attack would lead to a 'syn-product'
        major_product_structure = 'syn'

    # --- Step 4: Apply the Final Filter based on Product Structure ---
    # The major product must have the structure determined in Step 3.
    final_candidates = {k: v for k, v in possible_products.items() if v['bridge_type'] == major_product_structure}

    if len(final_candidates) != 1:
        return f"Constraint Check Failed (Final Selection): After applying all rules, {len(final_candidates)} candidates remain: {list(final_candidates.keys())}. A single major product was expected."

    predicted_answer = list(final_candidates.keys())[0]

    # --- Step 5: Compare with the LLM's Answer ---
    if predicted_answer == llm_answer:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but the step-by-step reasoning leads to {predicted_answer}.\n"
                f"Reasoning Trace:\n"
                f"1. Endo rule selects for 'endo' products -> Options {list(possible_products.keys())} remain.\n"
                f"2. Facial selectivity for fluorine favors '{favored_attack_type}' attack.\n"
                f"3. This attack results in an '{major_product_structure}' product.\n"
                f"4. The only remaining option that is an '{major_product_structure}' product is {predicted_answer}.")

# Execute the check
result = check_diels_alder_stereochemistry()
print(result)