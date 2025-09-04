def check_hplc_answer():
    """
    This function checks the correctness of the HPLC peak prediction by modeling the products
    and the separation principles of normal-phase and chiral HPLC.
    """

    # Step 1: Define all unique stereoisomers produced from both reactions.
    # Each product is represented as a dictionary containing its base name,
    # specific stereochemistry, and its relationship group for separation purposes.

    # Reaction I products: Two diastereomers
    # (S)-5-methoxyhexan-3-one -> (3R,5S)- and (3S,5S)-5-methoxyhexan-3-ol
    products_reaction_I = [
        {'name': '5-methoxyhexan-3-ol', 'stereochem': '(3R, 5S)', 'group_id': 'Rxn1_Diastereomer_A'},
        {'name': '5-methoxyhexan-3-ol', 'stereochem': '(3S, 5S)', 'group_id': 'Rxn1_Diastereomer_B'}
    ]

    # Reaction II products: One pair of enantiomers and one meso compound
    # pentane-2,4-dione -> (2R,4R)-, (2S,4S)-, and (2R,4S)-pentane-2,4-diol
    products_reaction_II = [
        # Enantiomeric pair (will co-elute in normal-phase HPLC)
        {'name': 'pentane-2,4-diol', 'stereochem': '(2R, 4R)', 'group_id': 'Rxn2_Enantiomeric_Pair'},
        {'name': 'pentane-2,4-diol', 'stereochem': '(2S, 4S)', 'group_id': 'Rxn2_Enantiomeric_Pair'},
        # Meso compound (diastereomer to the enantiomeric pair)
        {'name': 'pentane-2,4-diol', 'stereochem': '(2R, 4S)', 'group_id': 'Rxn2_Meso_Compound'}
    ]

    # Combine all products into a single list for analysis
    all_products = products_reaction_I + products_reaction_II

    # --- Analysis for Normal-Phase HPLC ---
    # Separates constitutional isomers and diastereomers. Enantiomers co-elute.
    # A peak's identity is defined by its constitutional isomer group ('name') and its diastereomer group ('group_id').
    normal_phase_peaks = set()
    for product in all_products:
        # Enantiomers share the same 'name' and 'group_id', so they will only be added to the set once.
        peak_signature = (product['name'], product['group_id'])
        normal_phase_peaks.add(peak_signature)
    
    calculated_normal_peaks = len(normal_phase_peaks)

    # --- Analysis for Chiral HPLC ---
    # Separates all unique stereoisomers (constitutional, diastereomers, and enantiomers).
    # A peak's identity is defined by its unique structure, represented here by 'name' and 'stereochem'.
    chiral_phase_peaks = set()
    for product in all_products:
        # Every unique stereoisomer has a unique combination of 'name' and 'stereochem'.
        peak_signature = (product['name'], product['stereochem'])
        chiral_phase_peaks.add(peak_signature)
        
    calculated_chiral_peaks = len(chiral_phase_peaks)

    # --- Verification ---
    # The answer 'D' states 5 peaks in chiral HPLC and 4 peaks in normal-phase HPLC.
    expected_chiral_peaks = 5
    expected_normal_peaks = 4

    if calculated_chiral_peaks == expected_chiral_peaks and calculated_normal_peaks == expected_normal_peaks:
        return "Correct"
    else:
        error_messages = []
        if calculated_chiral_peaks != expected_chiral_peaks:
            error_messages.append(
                f"Chiral HPLC peak count is incorrect. The answer states {expected_chiral_peaks}, but the analysis calculates {calculated_chiral_peaks}. "
                f"There are {len(all_products)} total unique stereoisomers, each of which should produce a peak in chiral HPLC."
            )
        if calculated_normal_peaks != expected_normal_peaks:
            error_messages.append(
                f"Normal-phase HPLC peak count is incorrect. The answer states {expected_normal_peaks}, but the analysis calculates {calculated_normal_peaks}. "
                f"The distinct groups are: {normal_phase_peaks}."
            )
        return "Incorrect. " + " ".join(error_messages)

# Execute the check and print the result
result = check_hplc_answer()
print(result)