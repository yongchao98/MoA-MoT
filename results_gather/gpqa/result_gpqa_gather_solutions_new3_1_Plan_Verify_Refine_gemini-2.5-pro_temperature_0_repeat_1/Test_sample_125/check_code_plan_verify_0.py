def check_hplc_peaks():
    """
    This function verifies the number of HPLC peaks based on the products of two chemical reactions.
    It follows the principles of stereochemistry and chromatography to calculate the expected number of peaks.
    """

    # --- Step 1: Analyze the products of Reaction I ---
    # Reactant: (S)-5-methoxyhexan-3-one. This is a chiral molecule with one stereocenter.
    # Reaction: Reduction of the ketone at C3 creates a new stereocenter.
    # Products: The original stereocenter (C5) is unchanged, while the new one (C3) can be R or S.
    # This results in two products: (3R, 5S)- and (3S, 5S)-5-methoxyhexan-3-ol.
    # Relationship: These two products are diastereomers.
    num_products_reaction_I = 2  # Two distinct diastereomers

    # --- Step 2: Analyze the products of Reaction II ---
    # Reactant: Pentane-2,4-dione. This is an achiral molecule.
    # Reaction: Reduction of both ketones creates two new stereocenters (C2, C4).
    # Products: The possible stereoisomers are (2R, 4R), (2S, 4S), and (2R, 4S).
    # Relationship:
    # - (2R, 4R) and (2S, 4S) are a pair of enantiomers.
    # - (2R, 4S) is a meso compound (achiral) and is a diastereomer of the enantiomeric pair.
    # This results in three distinct stereoisomers.
    num_enantiomeric_pairs_reaction_II = 1
    num_meso_compounds_reaction_II = 1

    # --- Step 3: Calculate peaks for Normal-Phase HPLC (Achiral Column) ---
    # A normal-phase (achiral) column separates diastereomers but NOT enantiomers.
    # It is assumed that products from Reaction I and Reaction II are structurally different
    # and will not co-elute.

    # Peaks from Reaction I:
    # The two products are diastereomers, so they have different physical properties and will be separated.
    normal_peaks_I = num_products_reaction_I  # 2 peaks

    # Peaks from Reaction II:
    # The enantiomeric pair will elute together as a single peak.
    # The meso compound is a diastereomer of the pair and will elute as a separate peak.
    normal_peaks_II = num_enantiomeric_pairs_reaction_II + num_meso_compounds_reaction_II  # 1 + 1 = 2 peaks

    # Total Normal-Phase Peaks:
    calculated_normal_peaks = normal_peaks_I + normal_peaks_II

    # --- Step 4: Calculate peaks for Chiral HPLC ---
    # A chiral HPLC column can separate all distinct stereoisomers, including enantiomers.

    # Peaks from Reaction I:
    # The two diastereomers will be separated.
    chiral_peaks_I = num_products_reaction_I  # 2 peaks

    # Peaks from Reaction II:
    # The two enantiomers in the pair will be resolved into two separate peaks.
    # The meso compound will give its own distinct peak.
    chiral_peaks_II = (num_enantiomeric_pairs_reaction_II * 2) + num_meso_compounds_reaction_II # (1*2) + 1 = 3 peaks

    # Total Chiral Peaks:
    calculated_chiral_peaks = chiral_peaks_I + chiral_peaks_II

    # --- Step 5: Verify the final answer ---
    # The question's answer is D: 5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC.
    expected_chiral_peaks = 5
    expected_normal_peaks = 4

    if calculated_chiral_peaks == expected_chiral_peaks and calculated_normal_peaks == expected_normal_peaks:
        return "Correct"
    else:
        reasoning = [
            "The provided answer is incorrect.",
            "Analysis:",
            f"1. Reaction I produces {num_products_reaction_I} diastereomers.",
            f"2. Reaction II produces 3 stereoisomers: 1 enantiomeric pair and 1 meso compound.",
            f"3. Normal-Phase HPLC (achiral) separates diastereomers but not enantiomers.",
            f"   - Peaks from Rxn I = {normal_peaks_I}",
            f"   - Peaks from Rxn II = {normal_peaks_II} (1 for enantiomeric pair, 1 for meso)",
            f"   - Total Calculated Normal-Phase Peaks = {calculated_normal_peaks}",
            f"4. Chiral HPLC separates all distinct stereoisomers.",
            f"   - Peaks from Rxn I = {chiral_peaks_I}",
            f"   - Peaks from Rxn II = {chiral_peaks_II} (2 for enantiomers, 1 for meso)",
            f"   - Total Calculated Chiral Peaks = {calculated_chiral_peaks}",
            "\nDiscrepancy:",
        ]
        if calculated_normal_peaks != expected_normal_peaks:
            reasoning.append(f"Normal-phase peaks do not match. Expected: {expected_normal_peaks}, Calculated: {calculated_normal_peaks}.")
        if calculated_chiral_peaks != expected_chiral_peaks:
            reasoning.append(f"Chiral peaks do not match. Expected: {expected_chiral_peaks}, Calculated: {calculated_chiral_peaks}.")
        
        return "\n".join(reasoning)

# Run the check
result = check_hplc_peaks()
print(result)