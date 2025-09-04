def check_organic_synthesis_answer():
    """
    Checks the correctness of the provided answer by applying the stereochemical rules
    of the three-step synthesis to the given options.
    """
    # The provided answer from the LLM
    llm_answer = 'D'

    # Verified properties of each isomer based on analysis of their SMILES strings and chemical structure.
    # 'diester': cis/trans relationship of the two ester groups.
    # 'precursor': Stereochemistry of the first Diels-Alder (syn/anti). 'anti' is major.
    # 'addition': Stereochemistry of the second Diels-Alder (endo/exo). 'exo' is major.
    isomer_properties = {
        'A': {'diester': 'trans', 'precursor': 'anti', 'addition': 'exo'},
        'B': {'diester': 'cis',   'precursor': 'anti', 'addition': 'endo'},
        'C': {'diester': 'cis',   'precursor': 'syn',  'addition': 'exo'},
        'D': {'diester': 'cis',   'precursor': 'anti', 'addition': 'exo'}
    }

    # --- Apply Chemical Rules to Find the Correct Major Product ---

    # Rule 1: Product must have a cis-diester configuration (from maleic anhydride).
    candidates = {key: props for key, props in isomer_properties.items() if props['diester'] == 'cis'}

    # Rule 2: Product must derive from the major 'anti' precursor (from 1st Diels-Alder).
    candidates = {key: props for key, props in candidates.items() if props['precursor'] == 'anti'}

    # Rule 3: Product must be the major 'exo' adduct (from 2nd Diels-Alder).
    candidates = {key: props for key, props in candidates.items() if props['addition'] == 'exo'}

    # --- Evaluate the Result ---

    # There should be exactly one isomer that satisfies all conditions for the major product.
    if len(candidates) != 1:
        return f"Checking Error: The chemical rules led to {len(candidates)} products instead of one."

    correct_answer = list(candidates.keys())[0]

    if llm_answer == correct_answer:
        return "Correct"
    else:
        # Determine why the LLM's answer is incorrect.
        incorrect_props = isomer_properties.get(llm_answer)
        if not incorrect_props:
            return f"Incorrect. The provided answer '{llm_answer}' is not a valid option."

        if incorrect_props['diester'] != 'cis':
            return (f"Incorrect. The provided answer '{llm_answer}' has a 'trans' diester configuration. "
                    "The reaction uses maleic anhydride, which requires the product to have a 'cis' configuration.")

        if incorrect_props['precursor'] != 'anti':
            return (f"Incorrect. The provided answer '{llm_answer}' derives from the minor 'syn' precursor of the first "
                    "Diels-Alder reaction. The major product forms from the sterically favored 'anti' precursor.")

        if incorrect_props['addition'] != 'exo':
            return (f"Incorrect. The provided answer '{llm_answer}' is the 'endo' adduct of the second "
                    "Diels-Alder reaction. The major product is the sterically favored 'exo' adduct.")
        
        return f"Incorrect. The correct answer is {correct_answer}, but the provided answer is {llm_answer}."

# Run the check and print the result.
result = check_organic_synthesis_answer()
print(result)