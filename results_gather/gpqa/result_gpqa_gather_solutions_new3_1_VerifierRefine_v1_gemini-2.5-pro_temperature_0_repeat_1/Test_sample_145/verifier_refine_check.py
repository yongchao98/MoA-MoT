def check_diels_alder_product():
    """
    Checks the correctness of the provided answer for the Diels-Alder reaction
    between 5-fluorocyclopenta-1,3-diene and maleic anhydride.
    """
    # The final answer from the LLM analysis to be checked.
    provided_answer = 'B'

    # Define the properties of the four options based on their IUPAC stereodescriptors.
    # Core Stereochemistry: 'endo' vs. 'exo'
    # Bridge Stereochemistry: 'syn' (F and anhydride on same side) vs. 'anti' (opposite sides)
    options = {
        'A': {'name': '(3aR,4R,7S,7aS,8r)-...', 'core_stereochem': 'exo', 'bridge_stereochem': 'anti'},
        'B': {'name': '(3aR,4S,7R,7aS,8r)-...', 'core_stereochem': 'endo', 'bridge_stereochem': 'anti'},
        'C': {'name': '(3aR,4R,7S,7aS,8s)-...', 'core_stereochem': 'exo', 'bridge_stereochem': 'syn'},
        'D': {'name': '(3aR,4S,7R,7aS,8s)-...', 'core_stereochem': 'endo', 'bridge_stereochem': 'syn'},
    }

    # --- Step 1: Apply the Endo Rule ---
    # The major kinetic product is the 'endo' adduct.
    preferred_core = 'endo'
    
    # --- Step 2: Determine Facial Selectivity ---
    # For a C5-fluoro substituent, electronic effects favor 'syn'-facial attack.
    # 'syn'-attack leads to an 'anti'-adduct (F is anti to the anhydride ring).
    preferred_bridge = 'anti'

    # --- Step 3: Identify the Correct Product ---
    # The major product should be the 'endo, anti' isomer.
    correct_option_key = None
    for key, properties in options.items():
        if properties['core_stereochem'] == preferred_core and properties['bridge_stereochem'] == preferred_bridge:
            correct_option_key = key
            break
            
    if correct_option_key is None:
        return "Logic error: No option matches the predicted 'endo, anti' stereochemistry."

    # --- Step 4: Check the Provided Answer ---
    if provided_answer == correct_option_key:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{provided_answer}' is incorrect. "
            f"The analysis points to '{correct_option_key}' as the major product.\n\n"
            "Reasoning:\n"
            f"1. Endo Rule: The reaction favors the '{preferred_core}' adduct. This eliminates options A and C, leaving B and D.\n"
            "2. Facial Selectivity: Electronic effects of the fluorine substituent favor 'syn'-attack.\n"
            f"3. Product Stereochemistry: 'Syn'-attack results in an '{preferred_bridge}'-adduct, where the fluorine is on the opposite side of the anhydride ring. This corresponds to the '8r' descriptor.\n"
            f"4. Conclusion: The major product is the 'endo, anti' isomer, which is option {correct_option_key}."
        )
        return reason

# Execute the check
result = check_diels_alder_product()
print(result)