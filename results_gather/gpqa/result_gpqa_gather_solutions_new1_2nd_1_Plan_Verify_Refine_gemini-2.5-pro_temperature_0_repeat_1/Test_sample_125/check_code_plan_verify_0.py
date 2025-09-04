def check_hplc_peaks_answer():
    """
    This function checks the correctness of the answer to the organic chemistry HPLC problem.
    It programmatically determines the number of products and their stereochemical relationships
    for each reaction, and then calculates the expected number of peaks for both normal-phase
    and chiral HPLC.
    """

    # --- Step 1: Analyze Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one. It has one existing stereocenter (C5).
    # Reaction: Reduction of a ketone (C3) to an alcohol, creating a new stereocenter (C3).
    # The original stereocenter is unaffected. The new one can be R or S.
    # Products: (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol.
    # Relationship: These are diastereomers.
    # Conclusion: Reaction I produces 2 unique, diastereomeric compounds.
    products_from_reaction_1 = {
        "description": "Two diastereomers",
        "unique_compounds": 2
    }

    # --- Step 2: Analyze Reaction II ---
    # Reactant: Pentane-2,4-dione. It is achiral and symmetric.
    # Reaction: Reduction of both ketones (C2, C4) to alcohols, creating two new stereocenters.
    # Products: (2R, 4R)-pentane-2,4-diol, (2S, 4S)-pentane-2,4-diol, and (2R, 4S)-pentane-2,4-diol.
    # Relationship: (2R,4R) and (2S,4S) are an enantiomeric pair. (2R,4S) is a meso compound.
    # Conclusion: Reaction II produces 3 unique stereoisomers.
    products_from_reaction_2 = {
        "description": "One enantiomeric pair and one meso compound",
        "enantiomeric_pairs": 1,
        "meso_compounds": 1,
        "unique_compounds": 3  # (2 enantiomers + 1 meso)
    }

    # --- Step 3: Calculate Peaks for Normal-Phase HPLC (Achiral) ---
    # Normal-phase HPLC separates diastereomers but not enantiomers.
    # It is assumed that products from Reaction I and Reaction II are structurally different
    # and will not co-elute.
    
    # Peaks from Reaction I: The 2 diastereomers will be separated.
    normal_peaks_rxn1 = products_from_reaction_1["unique_compounds"]
    
    # Peaks from Reaction II: The enantiomeric pair co-elutes as 1 peak. The meso compound is a
    # diastereomer to the pair and gives a separate peak.
    normal_peaks_rxn2 = products_from_reaction_2["enantiomeric_pairs"] + products_from_reaction_2["meso_compounds"]
    
    calculated_normal_peaks = normal_peaks_rxn1 + normal_peaks_rxn2

    # --- Step 4: Calculate Peaks for Chiral HPLC ---
    # Chiral HPLC separates all unique stereoisomers, including enantiomers.
    # The total number of peaks is the total number of unique compounds.
    calculated_chiral_peaks = products_from_reaction_1["unique_compounds"] + products_from_reaction_2["unique_compounds"]

    # --- Step 5: Verify the Provided Answer ---
    # The provided answer is A, which corresponds to 5 chiral peaks and 4 normal peaks.
    provided_answer = {
        "chiral": 5,
        "normal": 4
    }

    # Check if the calculated values match the provided answer.
    if calculated_chiral_peaks == provided_answer["chiral"] and calculated_normal_peaks == provided_answer["normal"]:
        return "Correct"
    else:
        error_message = []
        if calculated_chiral_peaks != provided_answer["chiral"]:
            error_message.append(
                f"Incorrect number of chiral HPLC peaks. The answer states {provided_answer['chiral']}, but the calculation shows {calculated_chiral_peaks}. "
                f"Reason: There are {products_from_reaction_1['unique_compounds']} products from Rxn I and {products_from_reaction_2['unique_compounds']} from Rxn II, totaling {calculated_chiral_peaks} unique stereoisomers, all of which are separated by a chiral column."
            )
        if calculated_normal_peaks != provided_answer["normal"]:
            error_message.append(
                f"Incorrect number of normal-phase HPLC peaks. The answer states {provided_answer['normal']}, but the calculation shows {calculated_normal_peaks}. "
                f"Reason: Rxn I gives {normal_peaks_rxn1} peaks (diastereomers). Rxn II gives {normal_peaks_rxn2} peaks (enantiomeric pair co-elutes, meso is separate). Total is {calculated_normal_peaks}."
            )
        return "\n".join(error_message)

# Execute the check and print the result
result = check_hplc_peaks_answer()
print(result)