def check_hplc_answer():
    """
    This function checks the correctness of the provided answer by analyzing the stereochemical outcomes
    of the reactions and applying the principles of normal-phase and chiral HPLC.
    """

    # --- Step 1: Analyze the products of each reaction ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one
    # The starting material has one chiral center at C5 (S-configuration).
    # The reduction of the ketone at C3 creates a new chiral center.
    # Since the starting material is enantiomerically pure, the two products formed,
    # (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol, are diastereomers.
    # Diastereomers are different compounds with different physical properties.
    num_stereoisomers_rxn1 = 2

    # Reaction II: Reduction of pentane-2,4-dione
    # The starting material is achiral.
    # The reduction of the two ketones at C2 and C4 creates two new chiral centers.
    # The possible products are (2R, 4R), (2S, 4S), and a meso compound (2R, 4S).
    # (2R, 4R) and (2S, 4S) form an enantiomeric pair.
    # The meso compound is a diastereomer of the enantiomeric pair.
    # This results in a total of 3 unique stereoisomers.
    num_stereoisomers_rxn2 = 3

    # --- Step 2: Calculate the expected number of peaks for each HPLC type ---

    # Chiral HPLC Analysis:
    # A chiral HPLC column can separate all unique stereoisomers, including enantiomers.
    # The total number of peaks will be the sum of all unique stereoisomers produced.
    # We assume no accidental co-elution between products from the two different reactions.
    expected_chiral_peaks = num_stereoisomers_rxn1 + num_stereoisomers_rxn2

    # Normal-phase HPLC Analysis:
    # A normal-phase (achiral) HPLC column separates compounds based on differences in
    # structure and polarity, such as constitutional isomers and diastereomers.
    # It cannot separate enantiomers, which elute together as a single peak.
    
    # Peaks from Reaction I: The two products are diastereomers, so they will be separated.
    peaks_from_rxn1 = 2
    
    # Peaks from Reaction II: The enantiomeric pair will elute as one peak. The meso compound,
    # being a diastereomer of the pair, will elute as a separate peak.
    peaks_from_rxn2 = 2 # (1 for the enantiomeric pair + 1 for the meso compound)

    # Total peaks: We assume the products of Reaction I (methoxyhexanols) are structurally
    # different enough from the products of Reaction II (pentanediols) to be separated,
    # which is a standard assumption for such problems.
    expected_normal_phase_peaks = peaks_from_rxn1 + peaks_from_rxn2

    # --- Step 3: Compare the calculated results with the provided answer (Option C) ---
    
    answer_chiral_peaks = 5
    answer_normal_phase_peaks = 4

    if expected_chiral_peaks != answer_chiral_peaks:
        return (f"Incorrect: The number of chiral HPLC peaks is wrong. "
                f"Calculation: Reaction I produces 2 diastereomers. Reaction II produces 3 stereoisomers (an enantiomeric pair and a meso compound). "
                f"Total unique stereoisomers = 2 + 3 = 5. A chiral column should show 5 peaks. "
                f"The answer expects {answer_chiral_peaks} peaks, but the calculation gives {expected_chiral_peaks}.")

    if expected_normal_phase_peaks != answer_normal_phase_peaks:
        return (f"Incorrect: The number of normal-phase HPLC peaks is wrong. "
                f"Calculation: Reaction I's 2 diastereomers give 2 peaks. Reaction II's enantiomeric pair gives 1 peak and its meso compound gives 1 peak (total 2). "
                f"Total normal-phase peaks = 2 + 2 = 4. "
                f"The answer expects {answer_normal_phase_peaks} peaks, but the calculation gives {expected_normal_phase_peaks}.")

    return "Correct"

# Run the check
result = check_hplc_answer()
print(result)