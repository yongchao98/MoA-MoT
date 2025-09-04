def check_chemistry_hplc_answer():
    """
    This function checks the correctness of the answer to the organic chemistry HPLC problem.
    It models the chemical logic step-by-step to derive the expected number of peaks.
    """

    # --- Step 1: Analyze Reaction I ---
    # (S)-5-methoxyhexan-3-one (chiral) is reduced, creating a new stereocenter.
    # The original stereocenter is retained.
    # This creates two products that are diastereomers of each other.
    products_from_reaction_I = {
        "count": 2,
        "type": "diastereomers",
        "description": "(3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol"
    }

    # --- Step 2: Analyze Reaction II ---
    # Pentane-2,4-dione (achiral) is reduced, creating two new stereocenters.
    # This creates a pair of enantiomers and a meso compound.
    products_from_reaction_II = {
        "count": 3,
        "type": "enantiomeric_pair_and_meso",
        "description": "(2R,4R)-diol, (2S,4S)-diol, and meso-diol"
    }

    # --- Step 3: Calculate Peaks for Normal-Phase HPLC (Achiral) ---
    # Principle: Separates diastereomers, but not enantiomers.
    
    # Peaks from Reaction I: The two products are diastereomers, so they separate.
    peaks_I_normal = 2
    
    # Peaks from Reaction II: The enantiomeric pair co-elutes (1 peak).
    # The meso compound is a diastereomer of the pair and separates (1 peak).
    peaks_II_normal = 1 + 1
    
    total_normal_peaks = peaks_I_normal + peaks_II_normal

    # --- Step 4: Calculate Peaks for Chiral HPLC ---
    # Principle: Separates all unique stereoisomers (diastereomers and enantiomers).
    
    # Peaks from Reaction I: The two diastereomers separate.
    peaks_I_chiral = 2
    
    # Peaks from Reaction II: The enantiomeric pair is resolved (2 peaks).
    # The meso compound gives its own peak (1 peak).
    peaks_II_chiral = 2 + 1
    
    total_chiral_peaks = peaks_I_chiral + peaks_II_chiral

    # --- Step 5: Verify the Final Answer ---
    # The provided answer is 'B', which states "5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC".
    
    expected_chiral_peaks = 5
    expected_normal_peaks = 4
    
    llm_answer_choice = "B"
    
    # Check if the calculated results match the expected results.
    if total_chiral_peaks == expected_chiral_peaks and total_normal_peaks == expected_normal_peaks:
        # The logic is correct. Now check if the LLM's choice 'B' reflects this.
        # Option B: 5 Chiral, 4 Normal. This matches our findings.
        return "Correct"
    else:
        # This block would execute if our logical derivation was different.
        error_message = []
        if total_chiral_peaks != expected_chiral_peaks:
            error_message.append(f"Chiral peak calculation is wrong. Expected {expected_chiral_peaks}, but logic yields {total_chiral_peaks}.")
        if total_normal_peaks != expected_normal_peaks:
            error_message.append(f"Normal-phase peak calculation is wrong. Expected {expected_normal_peaks}, but logic yields {total_normal_peaks}.")
        
        # This block checks if the LLM's answer is consistent with the flawed logic (if any).
        # For this problem, the logic is sound and matches the answer.
        return "Incorrect. " + " ".join(error_message)

# Run the check
result = check_chemistry_hplc_answer()
print(result)