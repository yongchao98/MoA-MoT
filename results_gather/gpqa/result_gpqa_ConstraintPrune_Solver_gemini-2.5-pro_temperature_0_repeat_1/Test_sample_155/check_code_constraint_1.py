def check_answer():
    """
    This function checks the correctness of the answer by simulating the stereochemical
    outcomes and the principles of achiral vs. chiral chromatography.
    """

    # --- Define Stereoisomers ---
    # We represent stereoisomers with a simple dictionary.
    # The 'id' helps define stereochemical relationships:
    # - Enantiomers have opposite non-zero ids (e.g., 1 and -1).
    # - A meso compound has an id of 0.
    # - Diastereomers have different absolute ids.
    meso_diol = {'name': '(4R,5S)-octane-4,5-diol', 'id': 0}
    enantiomer_RR = {'name': '(4R,5R)-octane-4,5-diol', 'id': 1}
    enantiomer_SS = {'name': '(4S,5S)-octane-4,5-diol', 'id': -1}

    # --- Simulate Reactions ---
    # Reaction 1: Anti-dihydroxylation of a trans-alkene yields a meso compound.
    products_reaction_1 = [meso_diol]

    # Reaction 2: Anti-dihydroxylation of a cis-alkene yields a racemic mixture.
    products_reaction_2 = [enantiomer_RR, enantiomer_SS]

    # --- Combine Products ---
    combined_mixture = products_reaction_1 + products_reaction_2
    
    # --- Simulate HPLC Analysis ---

    # 1. Standard (achiral) HPLC
    # In achiral HPLC, enantiomers co-elute (are not resolved). Diastereomers are resolved.
    # We can simulate this by grouping compounds that are enantiomers.
    # A simple way is to use the absolute value of the 'id' as a key for a peak.
    # The meso compound (id=0) is one peak. The enantiomers (id=1, id=-1) form another peak.
    standard_hplc_peaks = set()
    for compound in combined_mixture:
        peak_identifier = abs(compound['id'])
        standard_hplc_peaks.add(peak_identifier)
    num_standard_peaks = len(standard_hplc_peaks)

    # 2. Chiral HPLC
    # In chiral HPLC, all unique stereoisomers are resolved.
    # The number of peaks is simply the number of unique compounds in the mixture.
    num_chiral_peaks = len(combined_mixture)

    # --- Verify the Answer ---
    # The proposed answer (A) states 2 peaks in standard HPLC and 3 peaks in chiral HPLC.
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    if num_standard_peaks != expected_standard_peaks:
        return (f"Incorrect. The number of peaks in standard HPLC is wrong. "
                f"Expected: {expected_standard_peaks}, but the analysis shows there should be {num_standard_peaks}. "
                f"Reason: The mixture contains one meso compound and a pair of enantiomers. "
                f"Standard HPLC separates diastereomers (the meso compound from the enantiomeric pair) but does not resolve the enantiomers, resulting in 2 peaks.")

    if num_chiral_peaks != expected_chiral_peaks:
        return (f"Incorrect. The number of peaks in chiral HPLC is wrong. "
                f"Expected: {expected_chiral_peaks}, but the analysis shows there should be {num_chiral_peaks}. "
                f"Reason: The mixture contains three unique stereoisomers (one meso, two enantiomers). "
                f"Chiral HPLC resolves all of them, resulting in 3 peaks.")

    return "Correct"

# Run the check
result = check_answer()
print(result)