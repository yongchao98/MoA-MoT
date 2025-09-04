def check_hplc_peaks_answer():
    """
    This function verifies the number of HPLC peaks from the combined products of two reactions.
    It models the products and the separation principles of normal-phase and chiral HPLC.
    """

    # --- Step 1: Define the products from both reactions ---

    # Reaction I: (S)-5-methoxyhexan-3-one -> (3R,5S)- and (3S,5S)-5-methoxyhexan-3-ol
    # These two products are diastereomers. They are not enantiomers of each other.
    # They will be separated by both normal and chiral HPLC.
    # We assign them unique 'enantiomer_group' IDs to reflect they are not an enantiomeric pair.
    products_reaction_1 = [
        {
            'name': '(3R, 5S)-5-methoxyhexan-3-ol',
            'base_structure': '5-methoxyhexan-3-ol',
            'enantiomer_group': 'Rxn1_Diastereomer_1'
        },
        {
            'name': '(3S, 5S)-5-methoxyhexan-3-ol',
            'base_structure': '5-methoxyhexan-3-ol',
            'enantiomer_group': 'Rxn1_Diastereomer_2'
        }
    ]

    # Reaction II: Pentane-2,4-dione -> (2R,4R)-, (2S,4S)-, and meso-pentane-2,4-diol
    # This reaction produces a pair of enantiomers and a meso compound.
    # The enantiomers get the same 'enantiomer_group' ID because they co-elute on a normal-phase column.
    # The meso compound is a diastereomer to the enantiomers and gets its own unique ID.
    products_reaction_2 = [
        {
            'name': '(2R, 4R)-pentane-2,4-diol',
            'base_structure': 'pentane-2,4-diol',
            'enantiomer_group': 'Rxn2_Enantiomeric_Pair'
        },
        {
            'name': '(2S, 4S)-pentane-2,4-diol',
            'base_structure': 'pentane-2,4-diol',
            'enantiomer_group': 'Rxn2_Enantiomeric_Pair'
        },
        {
            'name': 'meso-pentane-2,4-diol',
            'base_structure': 'pentane-2,4-diol',
            'enantiomer_group': 'Rxn2_Meso_Compound'
        }
    ]

    # The final mixture contains all products
    all_products = products_reaction_1 + products_reaction_2

    # --- Step 2: Calculate expected peaks for each HPLC type ---

    # Chiral HPLC separates all unique stereoisomers.
    # The number of peaks is the total number of unique product molecules.
    num_chiral_peaks = len(all_products)

    # Normal-Phase HPLC separates constitutional isomers and diastereomers, but not enantiomers.
    # We count peaks by identifying unique combinations of 'base_structure' and 'enantiomer_group'.
    # Using a set automatically handles uniqueness.
    normal_phase_peak_identifiers = set()
    for product in all_products:
        identifier = (product['base_structure'], product['enantiomer_group'])
        normal_phase_peak_identifiers.add(identifier)
    
    num_normal_peaks = len(normal_phase_peak_identifiers)

    # --- Step 3: Compare calculated results with the LLM's answer ---
    
    # The LLM's answer states: 5 peaks in chiral HPLC and 4 peaks in normal-phase HPLC.
    llm_chiral_peaks_answer = 5
    llm_normal_peaks_answer = 4

    if num_chiral_peaks == llm_chiral_peaks_answer and num_normal_peaks == llm_normal_peaks_answer:
        return "Correct"
    else:
        error_message = "The answer is incorrect.\n"
        if num_chiral_peaks != llm_chiral_peaks_answer:
            error_message += (f"Constraint failure (Chiral HPLC): The number of peaks should be {num_chiral_peaks}, "
                              f"but the answer claims {llm_chiral_peaks_answer}. Chiral HPLC separates all {len(all_products)} "
                              f"unique stereoisomers.\n")
        if num_normal_peaks != llm_normal_peaks_answer:
            error_message += (f"Constraint failure (Normal-Phase HPLC): The number of peaks should be {num_normal_peaks}, "
                              f"but the answer claims {llm_normal_peaks_answer}. Normal-phase HPLC separates the 2 diastereomers from Rxn I (2 peaks) "
                              f"and separates the meso compound from the enantiomeric pair in Rxn II (2 peaks), for a total of 4 peaks.\n")
        return error_message

# Run the check and print the result.
result = check_hplc_peaks_answer()
print(result)