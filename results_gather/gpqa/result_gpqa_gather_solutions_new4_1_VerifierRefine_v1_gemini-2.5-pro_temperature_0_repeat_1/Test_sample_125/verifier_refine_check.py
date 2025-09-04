def check_hplc_answer():
    """
    This function programmatically verifies the number of HPLC peaks based on the products of two chemical reactions.
    It follows the principles of stereochemistry and chromatography to calculate the expected number of peaks
    and compares them to the provided answer.
    """

    # --- Step 1: Analyze the products of Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one is a chiral molecule with one stereocenter (C5).
    # Reaction: Reduction of the ketone at C3 creates a new stereocenter.
    # Products: Since the starting material is chiral, the reduction creates two diastereomers:
    # (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol.
    # Conclusion: Reaction I produces 2 distinct stereoisomers (a pair of diastereomers).
    products_from_reaction_1 = {
        "count": 2,
        "type": "diastereomers"
    }

    # --- Step 2: Analyze the products of Reaction II ---
    # Reactant: Pentane-2,4-dione is an achiral molecule.
    # Reaction: Reduction of both ketones at C2 and C4 creates two new stereocenters.
    # Products: Since the starting material is achiral, all possible stereoisomers are formed.
    # 1. (2R, 4R)-pentane-2,4-diol and (2S, 4S)-pentane-2,4-diol (a pair of enantiomers).
    # 2. (2R, 4S)-pentane-2,4-diol (a meso compound, which is achiral).
    # Conclusion: Reaction II produces 3 distinct stereoisomers.
    products_from_reaction_2 = {
        "count": 3,
        "enantiomeric_pairs": 1,  # The (2R,4R) and (2S,4S) pair
        "meso_compounds": 1       # The (2R,4S) compound
    }

    # --- Step 3: Calculate peaks for Normal-Phase HPLC (achiral column) ---
    # Principle: Separates diastereomers and constitutional isomers, but NOT enantiomers.
    
    # Peaks from Reaction I: The 2 diastereomers have different physical properties and will be separated.
    normal_peaks_rxn1 = 2
    
    # Peaks from Reaction II:
    # The enantiomeric pair co-elutes as a single peak.
    # The meso compound is a diastereomer of the enantiomeric pair and will be separated.
    normal_peaks_rxn2 = products_from_reaction_2["enantiomeric_pairs"] + products_from_reaction_2["meso_compounds"]  # 1 + 1 = 2
    
    # Total normal-phase peaks. We assume products from Rxn I and II are structurally different and will not co-elute.
    calculated_normal_peaks = normal_peaks_rxn1 + normal_peaks_rxn2

    # --- Step 4: Calculate peaks for Chiral HPLC ---
    # Principle: Separates all unique stereoisomers, including both diastereomers and enantiomers.
    
    # Peaks from Reaction I: The 2 diastereomers will be separated.
    chiral_peaks_rxn1 = products_from_reaction_1["count"]
    
    # Peaks from Reaction II:
    # The enantiomeric pair is resolved into two separate peaks.
    # The meso compound gives its own distinct peak.
    chiral_peaks_rxn2 = (products_from_reaction_2["enantiomeric_pairs"] * 2) + products_from_reaction_2["meso_compounds"] # (1 * 2) + 1 = 3
    
    # Total chiral peaks.
    calculated_chiral_peaks = chiral_peaks_rxn1 + chiral_peaks_rxn2

    # --- Step 5: Verify the provided answer ---
    # The final answer provided in the prompt is <<<D>>>.
    # Let's define the options from the question.
    options = {
        'A': {'chiral': 3, 'normal': 2},
        'B': {'chiral': 3, 'normal': 3},
        'C': {'chiral': 4, 'normal': 2},
        'D': {'chiral': 5, 'normal': 4}
    }
    
    final_answer_choice = 'D'
    answer_values = options[final_answer_choice]
    
    # Compare calculated values with the values from the chosen answer option.
    is_correct = True
    error_messages = []

    if calculated_chiral_peaks != answer_values['chiral']:
        is_correct = False
        error_messages.append(
            f"Incorrect number of chiral HPLC peaks. "
            f"The code calculated {calculated_chiral_peaks} peaks, but the answer D states {answer_values['chiral']} peaks. "
            f"Calculation: Reaction I (2 diastereomers) -> 2 peaks. Reaction II (1 enantiomeric pair + 1 meso) -> 3 peaks. Total = 2 + 3 = 5."
        )

    if calculated_normal_peaks != answer_values['normal']:
        is_correct = False
        error_messages.append(
            f"Incorrect number of normal-phase HPLC peaks. "
            f"The code calculated {calculated_normal_peaks} peaks, but the answer D states {answer_values['normal']} peaks. "
            f"Calculation: Reaction I (2 diastereomers) -> 2 peaks. Reaction II (1 enantiomeric pair + 1 meso) -> 2 peaks. Total = 2 + 2 = 4."
        )

    if is_correct:
        return "Correct"
    else:
        return "\n".join(error_messages)

# Execute the check and print the result
result = check_hplc_answer()
print(result)