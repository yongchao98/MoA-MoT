def check_hplc_peaks():
    """
    This function checks the correctness of the LLM's answer by modeling the chemical reactions and chromatographic separations.

    It determines the number of stereoisomers produced in each reaction and then calculates the number of expected peaks
    in both normal-phase and chiral-phase HPLC based on the principles of stereoisomer separation.
    """

    # --- Step 1: Define the products of each reaction based on stereochemistry principles ---

    # Reaction I: (S)-5-methoxyhexan-3-one (chiral) is reduced by LAH (achiral).
    # This creates a new stereocenter at C3. The original C5 stereocenter is unaffected.
    # The products are (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol.
    # These two products are diastereomers.
    # We represent them as unique strings for counting purposes.
    products_reaction_I = [
        "diastereomer_1_from_rxn_I",  # e.g., (3R, 5S)-5-methoxyhexan-3-ol
        "diastereomer_2_from_rxn_I"   # e.g., (3S, 5S)-5-methoxyhexan-3-ol
    ]
    num_products_I = len(products_reaction_I)

    # Reaction II: Pentane-2,4-dione (achiral) is reduced by NaBH4 (achiral).
    # This creates two new stereocenters at C2 and C4.
    # The possible products are (2R, 4R), (2S, 4S), and (2R, 4S).
    # (2R, 4R) and (2S, 4S) are a pair of enantiomers.
    # (2R, 4S) is an achiral meso compound.
    # This gives a total of 3 unique stereoisomers.
    products_reaction_II = {
        "enantiomeric_pair": ["(2R, 4R)-pentane-2,4-diol", "(2S, 4S)-pentane-2,4-diol"],
        "meso_compound": ["(2R, 4S)-pentane-2,4-diol"]
    }
    num_products_II = len(products_reaction_II["enantiomeric_pair"]) + len(products_reaction_II["meso_compound"])

    # Total unique stereoisomers produced from both reactions
    total_unique_compounds = num_products_I + num_products_II

    # --- Step 2: Calculate the expected number of peaks for each HPLC type ---

    # Chiral HPLC: Separates all unique stereoisomers (both diastereomers and enantiomers).
    # The number of peaks is equal to the total number of unique compounds.
    # We assume no accidental co-elution between products of Reaction I and Reaction II,
    # which is a safe assumption given their structural differences (C7 ether-alcohol vs C5 diol).
    expected_chiral_peaks = total_unique_compounds

    # Normal-Phase HPLC: Separates compounds based on polarity.
    # - Diastereomers have different physical properties and will separate into different peaks.
    # - Enantiomers have identical physical properties in an achiral environment and will co-elute as a single peak.
    
    # Peaks from Reaction I: The two products are diastereomers, so they will separate.
    peaks_from_I_normal = num_products_I  # 2 peaks

    # Peaks from Reaction II: The enantiomeric pair co-elutes as 1 peak. The meso compound is a diastereomer
    # to the enantiomeric pair and will have a different polarity, giving a separate peak.
    peaks_from_II_normal = 1  # for the enantiomeric pair
    peaks_from_II_normal += 1 # for the meso compound
    
    # Total peaks in normal-phase HPLC
    expected_normal_phase_peaks = peaks_from_I_normal + peaks_from_II_normal

    # --- Step 3: Compare the calculated results with the LLM's answer ---
    
    # The LLM's answer corresponds to option A.
    llm_answer_chiral_peaks = 5
    llm_answer_normal_phase_peaks = 4

    errors = []
    if expected_chiral_peaks != llm_answer_chiral_peaks:
        errors.append(
            f"Incorrect number of peaks for Chiral HPLC. "
            f"Calculation: Reaction I produces {num_products_I} diastereomers. Reaction II produces {num_products_II} stereoisomers (1 enantiomeric pair + 1 meso compound). "
            f"Total unique compounds = {num_products_I} + {num_products_II} = {total_unique_compounds}. "
            f"Chiral HPLC separates all unique compounds, so it should show {expected_chiral_peaks} peaks. "
            f"The answer states {llm_answer_chiral_peaks} peaks."
        )

    if expected_normal_phase_peaks != llm_answer_normal_phase_peaks:
        errors.append(
            f"Incorrect number of peaks for Normal-Phase HPLC. "
            f"Calculation: Reaction I's {num_products_I} diastereomers give {peaks_from_I_normal} peaks. "
            f"Reaction II's enantiomeric pair gives 1 peak and its meso compound gives 1 peak, for a total of {peaks_from_II_normal} peaks. "
            f"Total peaks = {peaks_from_I_normal} + {peaks_from_II_normal} = {expected_normal_phase_peaks}. "
            f"The answer states {llm_answer_normal_phase_peaks} peaks."
        )

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check
result = check_hplc_peaks()
print(result)