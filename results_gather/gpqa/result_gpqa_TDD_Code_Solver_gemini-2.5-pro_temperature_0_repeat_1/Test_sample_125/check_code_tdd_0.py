def check_hplc_answer():
    """
    This function verifies the logic for determining the number of HPLC peaks
    from the products of two organic reactions.

    It independently calculates the expected number of peaks based on stereochemical principles
    and compares it to the provided answer's claim.
    """

    # --- Step 1: Analyze Reaction I ---
    # Starting material: (S)-5-methoxyhexan-3-one. It has one chiral center (C5).
    # Reaction: Reduction of the ketone at C3 to an alcohol. This creates a new chiral center at C3.
    # Since a new chiral center is formed in a molecule that already has one,
    # the products are two diastereomers: (3R, 5S) and (3S, 5S)-5-methoxyhexan-3-ol.
    # Diastereomers are chemically and physically distinct compounds.
    products_reaction_1 = 2  # Number of diastereomers

    # --- Step 2: Analyze Reaction II ---
    # Starting material: Pentane-2,4-dione. This is an achiral molecule (has a plane of symmetry).
    # Reaction: Reduction of both ketones to alcohols, forming pentane-2,4-diol.
    # This creates two new chiral centers (C2 and C4).
    # The possible stereoisomers are (2R,4R), (2S,4S), and (2R,4S).
    # (2R,4R) and (2S,4S) are a pair of enantiomers.
    # (2R,4S) is a meso compound (achiral, with a plane of symmetry).
    # Total unique stereoisomers from Reaction II = 3 (the two enantiomers + the meso compound).
    enantiomeric_pairs_reaction_2 = 1
    meso_compounds_reaction_2 = 1
    total_stereoisomers_reaction_2 = 2 * enantiomeric_pairs_reaction_2 + meso_compounds_reaction_2

    # --- Step 3: Analyze Chiral HPLC Separation ---
    # A chiral HPLC column can separate all unique stereoisomers (enantiomers and diastereomers).
    # The products from Reaction I and Reaction II are constitutionally different (C7 vs C5 compounds),
    # so they will definitely separate from each other.
    # Total peaks = (Products from Rxn I) + (Total stereoisomers from Rxn II)
    expected_chiral_peaks = products_reaction_1 + total_stereoisomers_reaction_2
    
    # --- Step 4: Analyze Normal-Phase HPLC Separation ---
    # A normal-phase (achiral) HPLC column separates constitutional isomers and diastereomers,
    # but it CANNOT separate enantiomers.
    # Peaks from Rxn I: The two products are diastereomers, so they will separate.
    peaks_from_rxn1_normal = products_reaction_1
    # Peaks from Rxn II: The enantiomeric pair will co-elute as one peak. The meso compound is a
    # diastereomer of the enantiomeric pair, so it will separate, giving a second peak.
    peaks_from_rxn2_normal = enantiomeric_pairs_reaction_2 + meso_compounds_reaction_2
    # Total peaks = (Peaks from Rxn I) + (Peaks from Rxn II)
    expected_normal_peaks = peaks_from_rxn1_normal + peaks_from_rxn2_normal

    # --- Step 5: Verify the Answer ---
    # The provided answer (Option A) claims 5 peaks in chiral HPLC and 4 peaks in normal-phase HPLC.
    answer_chiral_peaks = 5
    answer_normal_peaks = 4

    if expected_chiral_peaks != answer_chiral_peaks:
        return (f"Incorrect. The number of chiral HPLC peaks is wrong. "
                f"Calculation: {products_reaction_1} (from Rxn I) + {total_stereoisomers_reaction_2} (from Rxn II) = {expected_chiral_peaks}. "
                f"The answer claims {answer_chiral_peaks}.")

    if expected_normal_peaks != answer_normal_peaks:
        return (f"Incorrect. The number of normal-phase HPLC peaks is wrong. "
                f"Calculation: {peaks_from_rxn1_normal} (diastereomers from Rxn I) + {peaks_from_rxn2_normal} (enantiomeric pair + meso from Rxn II) = {expected_normal_peaks}. "
                f"The answer claims {answer_normal_peaks}.")

    return "Correct"

# Run the check
result = check_hplc_answer()
print(result)