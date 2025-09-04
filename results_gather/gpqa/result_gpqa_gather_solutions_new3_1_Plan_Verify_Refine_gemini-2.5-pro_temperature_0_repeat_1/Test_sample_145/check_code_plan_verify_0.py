def check_diels_alder_product():
    """
    Checks the correctness of the final answer for the Diels-Alder reaction between
    5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The function codifies the key chemical principles and checks if the provided
    answer aligns with the predicted major product.
    """
    # The final answer provided by the LLM analysis.
    llm_answer = 'A'

    # Step 1: Define the properties of each possible product based on IUPAC nomenclature.
    # There is a strong consensus in the provided analyses for these mappings:
    # - Ring Fusion: (3aR,4S,7R,7aS...) corresponds to the 'endo' adduct.
    # - Substituent Position: '8r' corresponds to 'anti' (F is opposite to the anhydride ring),
    #   and '8s' corresponds to 'syn' (F is on the same side as the anhydride ring).
    options = {
        'A': {'ring_fusion': 'endo', 'substituent_position': 'anti'},
        'B': {'ring_fusion': 'exo',  'substituent_position': 'anti'},
        'C': {'ring_fusion': 'endo', 'substituent_position': 'syn'},
        'D': {'ring_fusion': 'exo',  'substituent_position': 'syn'}
    }

    # Step 2: Apply the first principle of stereoselectivity: Endo/Exo Selectivity.
    # The Alder-endo rule states that for kinetic control, the 'endo' product is favored
    # due to secondary orbital interactions.
    favored_ring_fusion = 'endo'
    
    # Step 3: Apply the second principle: Facial Selectivity (Syn/Anti Attack).
    # For small, electronegative substituents (like Fluorine) at the C5 position of
    # cyclopentadiene, electronic effects favor 'syn'-facial attack over 'anti'-attack.
    favored_attack_pathway = 'syn-attack'

    # Step 4: Map the attack pathway to the final product's relative stereochemistry.
    # - A 'syn'-attack (dienophile approaches from the same side as F) results in the
    #   fluorine being 'anti' to the anhydride ring in the final product.
    # - An 'anti'-attack would result in a 'syn' product.
    if favored_attack_pathway == 'syn-attack':
        favored_substituent_position = 'anti'
    else: # favored_attack_pathway == 'anti-attack'
        favored_substituent_position = 'syn'

    # Step 5: Determine the predicted correct option based on the chemical principles.
    predicted_correct_option = None
    for option, properties in options.items():
        if (properties['ring_fusion'] == favored_ring_fusion and
            properties['substituent_position'] == favored_substituent_position):
            predicted_correct_option = option
            break

    # Step 6: Compare the predicted correct option with the LLM's final answer.
    if llm_answer == predicted_correct_option:
        return "Correct"
    else:
        reason = f"The provided answer '{llm_answer}' is incorrect. "
        reason += f"The analysis based on established chemical principles predicts option '{predicted_correct_option}'.\n"
        reason += f"Reasoning:\n"
        reason += f"1. Endo/Exo Selectivity: The major product should be the '{favored_ring_fusion}' adduct. This eliminates options B and D.\n"
        reason += f"2. Facial Selectivity: For a C5-Fluorine substituent, electronic effects favor '{favored_attack_pathway}'.\n"
        reason += f"3. Product Stereochemistry: This pathway leads to a product where the fluorine is '{favored_substituent_position}' relative to the anhydride ring.\n"
        reason += f"The only option that is both '{favored_ring_fusion}' and '{favored_substituent_position}' is option '{predicted_correct_option}'."
        return reason

# Execute the check and print the result.
result = check_diels_alder_product()
print(result)