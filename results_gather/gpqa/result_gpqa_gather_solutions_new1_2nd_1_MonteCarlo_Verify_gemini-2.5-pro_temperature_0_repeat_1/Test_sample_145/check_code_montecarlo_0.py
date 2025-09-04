def check_diels_alder_product():
    """
    Checks the correctness of the final answer for the Diels-Alder reaction of
    5-fluorocyclopenta-1,3-diene with maleic anhydride.

    The function verifies the answer against established stereochemical principles:
    1. The Alder Endo Rule (kinetic control favors the 'endo' product).
    2. Facial selectivity for 5-fluorocyclopentadiene (electronic effects favor 'syn-attack').
    3. The resulting product structure ('syn-attack' leads to the 'anti-product').
    4. The IUPAC nomenclature for the product ('8r' for 'anti', '8s' for 'syn').
    """

    # The final answer provided by the LLM to be checked.
    llm_answer = 'A'

    # Define the four possible products from the question.
    options = {
        'A': "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'C': "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        'D': "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
    }

    # --- Step 1: Apply the Alder Endo Rule ---
    # Principle: The major product is 'endo'.
    # The 'endo' adduct has the core stereochemistry (3aR,4S,7R,7aS).
    # The 'exo' adduct has the core stereochemistry (3aR,4R,7S,7aS).
    
    endo_core_stereochem = "(3aR,4S,7R,7aS"
    
    # Filter for 'endo' products.
    endo_products = {key: name for key, name in options.items() if name.startswith(endo_core_stereochem)}

    if llm_answer not in endo_products:
        return (f"Incorrect. The provided answer '{llm_answer}' corresponds to an 'exo' product. "
                "The Alder Endo Rule states that the major kinetic product of a Diels-Alder reaction "
                "is the 'endo' adduct.")

    # --- Step 2: Apply Facial Selectivity and Determine Product Structure ---
    # Principle: For 5-fluorocyclopentadiene, electronic effects favor 'syn-facial attack'.
    # Principle: A 'syn-facial attack' in an 'endo' orientation results in the 'anti-product'
    # (where the fluorine and anhydride ring are on opposite sides of the bicyclic system).
    
    # Therefore, the major product must be the 'anti-product'.
    favored_product_structure = 'anti'

    # --- Step 3: Match Product Structure to IUPAC Nomenclature ---
    # Principle: In this bicyclic system, the 'anti' position for the C8 substituent is designated '8r'.
    # The 'syn' position is designated '8s'.
    
    # Determine the predicted correct option key.
    predicted_key = None
    for key in endo_products:
        if '8r' in options[key]: # '8r' corresponds to the 'anti' product
            if favored_product_structure == 'anti':
                predicted_key = key
        elif '8s' in options[key]: # '8s' corresponds to the 'syn' product
            if favored_product_structure == 'syn':
                predicted_key = key

    if predicted_key is None:
        # This case should not be reached with the current logic.
        return "Error in checker: Could not determine a predicted product based on the rules."

    # --- Step 4: Final Verification ---
    if llm_answer == predicted_key:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is '{llm_answer}', but the correct answer should be '{predicted_key}'.\n"
                f"Reasoning:\n"
                f"1. The major product must be 'endo', which narrows the choices to A and C.\n"
                f"2. Electronic effects favor 'syn-facial attack', which results in the 'anti-product'.\n"
                f"3. The 'anti-product' is designated with the '8r' stereodescriptor.\n"
                f"Therefore, the correct product is the 'endo, anti' adduct, which is option {predicted_key}.")

# Execute the check and print the result.
result = check_diels_alder_product()
print(result)