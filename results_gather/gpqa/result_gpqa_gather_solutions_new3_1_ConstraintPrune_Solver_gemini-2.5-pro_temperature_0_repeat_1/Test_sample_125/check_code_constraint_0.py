def check_chemistry_hplc_answer():
    """
    This function checks the correctness of the answer by modeling the chemical reactions
    and the principles of HPLC separation.
    """

    # --- Step 1: Model the products of Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one (chiral, 1 stereocenter)
    # Reaction: Reduction of ketone at C3 creates a new stereocenter.
    # Outcome: The original stereocenter (C5) is unchanged. The new one (C3) can be R or S.
    # This creates two products: (3R, 5S) and (3S, 5S). These are diastereomers.
    products_reaction_1 = {
        "compound_A": {"name": "(3R, 5S)-5-methoxyhexan-3-ol", "group": "Rxn1_Diastereomer_1"},
        "compound_B": {"name": "(3S, 5S)-5-methoxyhexan-3-ol", "group": "Rxn1_Diastereomer_2"}
    }
    num_products_rxn1 = len(products_reaction_1)

    # --- Step 2: Model the products of Reaction II ---
    # Reactant: Pentane-2,4-dione (achiral)
    # Reaction: Reduction of both ketones creates two new stereocenters (C2, C4).
    # Outcome: Forms all possible stereoisomers.
    # (2R, 4R) and (2S, 4S) are an enantiomeric pair.
    # (2R, 4S) is a meso compound (diastereomer of the pair).
    products_reaction_2 = {
        "compound_C": {"name": "(2R, 4R)-pentane-2,4-diol", "group": "Rxn2_Enantiomeric_Pair"},
        "compound_D": {"name": "(2S, 4S)-pentane-2,4-diol", "group": "Rxn2_Enantiomeric_Pair"},
        "compound_E": {"name": "(2R, 4S)-pentane-2,4-diol", "group": "Rxn2_Meso_Compound"}
    }
    num_products_rxn2 = len(products_reaction_2)

    # --- Step 3: Combine all products ---
    all_products = {**products_reaction_1, **products_reaction_2}
    total_unique_stereoisomers = len(all_products)

    # --- Step 4: Calculate expected peaks for Chiral HPLC ---
    # Chiral HPLC separates all unique stereoisomers (enantiomers and diastereomers).
    # Therefore, the number of peaks is the total number of unique compounds.
    calculated_chiral_peaks = total_unique_stereoisomers

    # --- Step 5: Calculate expected peaks for Normal-Phase HPLC ---
    # Normal-phase HPLC separates diastereomers and constitutional isomers, but NOT enantiomers.
    # We can count the number of unique "groups" from our model.
    normal_phase_groups = set()
    for product_info in all_products.values():
        normal_phase_groups.add(product_info["group"])
    
    calculated_normal_peaks = len(normal_phase_groups)

    # --- Step 6: Get the proposed answer from the LLM ---
    # The final answer provided is <<<C>>>, which corresponds to:
    # C) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
    proposed_chiral_peaks = 5
    proposed_normal_peaks = 4

    # --- Step 7: Compare calculated results with the proposed answer ---
    error_messages = []
    if calculated_chiral_peaks != proposed_chiral_peaks:
        error_messages.append(
            f"Chiral HPLC peak count is incorrect. "
            f"Calculation shows {calculated_chiral_peaks} peaks (one for each of the {total_unique_stereoisomers} unique stereoisomers), "
            f"but the answer claims {proposed_chiral_peaks}."
        )
    
    if calculated_normal_peaks != proposed_normal_peaks:
        error_messages.append(
            f"Normal-Phase HPLC peak count is incorrect. "
            f"Calculation shows {calculated_normal_peaks} peaks (Rxn I: 2 diastereomers; Rxn II: 1 enantiomeric pair + 1 meso compound), "
            f"but the answer claims {proposed_normal_peaks}."
        )

    if not error_messages:
        return "Correct"
    else:
        return "Incorrect. " + " ".join(error_messages)

# Run the check
result = check_chemistry_hplc_answer()
print(result)