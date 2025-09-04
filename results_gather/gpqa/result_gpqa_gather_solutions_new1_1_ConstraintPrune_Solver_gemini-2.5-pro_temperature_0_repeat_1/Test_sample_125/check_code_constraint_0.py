def check_chemistry_hplc_problem():
    """
    This function checks the correctness of the answer to the organic chemistry HPLC problem.
    It codifies the analysis of the reaction products and the principles of chromatography.
    """

    # --- Step 1: Analyze the products of Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one (a single enantiomer, chiral)
    # Reaction: Reduction of a ketone creates a new stereocenter at C3.
    # The original stereocenter at C5 is unaffected.
    # Outcome: Formation of two products, (3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol.
    # Relationship: These are diastereomers.
    products_from_reaction_1 = {
        "description": "A pair of diastereomers",
        "count": 2
    }

    # --- Step 2: Analyze the products of Reaction II ---
    # Reactant: Pentane-2,4-dione (achiral)
    # Reaction: Reduction of two ketones creates two new stereocenters at C2 and C4.
    # Outcome: Formation of (2R, 4R)-, (2S, 4S)-, and (2R, 4S)-pentane-2,4-diol.
    # Relationship: (2R, 4R) and (2S, 4S) are an enantiomeric pair. (2R, 4S) is a meso compound.
    # In total, there are 3 unique stereoisomers.
    products_from_reaction_2 = {
        "description": "One enantiomeric pair and one meso compound",
        "enantiomeric_pairs": 1,
        "meso_compounds": 1,
        "total_stereoisomers": 3
    }

    # --- Step 3: Calculate the number of peaks in Normal-Phase HPLC ---
    # Principle: Separates diastereomers and constitutional isomers, but NOT enantiomers.
    # Peaks from Rxn I: The two products are diastereomers, so they will be separated.
    normal_hplc_peaks_rxn1 = products_from_reaction_1["count"]
    # Peaks from Rxn II: The enantiomeric pair co-elutes (1 peak). The meso compound is a
    # diastereomer of the pair and separates (1 peak).
    normal_hplc_peaks_rxn2 = products_from_reaction_2["enantiomeric_pairs"] + products_from_reaction_2["meso_compounds"]
    
    # Total peaks, assuming products from Rxn I and Rxn II are structurally different enough to separate,
    # which is a valid assumption (C7 ether-alcohol vs C5 diol).
    total_normal_hplc_peaks = normal_hplc_peaks_rxn1 + normal_hplc_peaks_rxn2

    # --- Step 4: Calculate the number of peaks in Chiral HPLC ---
    # Principle: Separates all unique stereoisomers, including enantiomers.
    # The total number of peaks is the total number of unique stereoisomers in the mixture.
    total_chiral_hplc_peaks = products_from_reaction_1["count"] + products_from_reaction_2["total_stereoisomers"]

    # --- Step 5: Compare with the provided answer ---
    # The question asks for the number of peaks in chiral and normal-phase HPLC.
    # The correct analysis yields:
    # Chiral HPLC: 5 peaks
    # Normal-phase HPLC: 4 peaks
    
    # The provided final answer is 'A', which corresponds to:
    # A) 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC
    
    expected_chiral_peaks = 5
    expected_normal_peaks = 4
    
    answer_chiral_peaks = 5
    answer_normal_peaks = 4

    if total_chiral_hplc_peaks == expected_chiral_peaks and total_normal_hplc_peaks == expected_normal_peaks:
        # The chemical analysis is correct. Now check if the provided answer matches.
        if answer_chiral_peaks == expected_chiral_peaks and answer_normal_peaks == expected_normal_peaks:
            return "Correct"
        else:
            return (f"The provided answer is incorrect. It states {answer_chiral_peaks} chiral peaks and "
                    f"{answer_normal_peaks} normal-phase peaks. The correct analysis shows there should be "
                    f"{expected_chiral_peaks} chiral peaks and {expected_normal_peaks} normal-phase peaks.")
    else:
        # This case would indicate an error in the logic of this checking script itself.
        return (f"Error in the checking script's logic. "
                f"Calculated Chiral Peaks: {total_chiral_hplc_peaks} (Expected: {expected_chiral_peaks}). "
                f"Calculated Normal Peaks: {total_normal_hplc_peaks} (Expected: {expected_normal_peaks}).")

# Execute the check
result = check_chemistry_hplc_problem()
print(result)