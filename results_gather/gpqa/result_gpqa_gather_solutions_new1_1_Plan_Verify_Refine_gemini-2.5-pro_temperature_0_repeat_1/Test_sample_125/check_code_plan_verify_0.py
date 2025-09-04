def check_hplc_peaks():
    """
    This function verifies the number of HPLC peaks based on the products of the two reactions.
    """
    # --- Step 1: Define the products and their stereochemical relationships ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one
    # The starting material is chiral. The reaction creates a new stereocenter.
    # The two products are diastereomers of each other.
    # Let's represent them as two distinct, non-enantiomeric compounds.
    rxn1_products = [
        {"id": "1A", "name": "(3R, 5S)-5-methoxyhexan-3-ol", "enantiomer_id": None},
        {"id": "1B", "name": "(3S, 5S)-5-methoxyhexan-3-ol", "enantiomer_id": None}
    ]

    # Reaction II: Reduction of pentane-2,4-dione
    # The starting material is achiral. The reaction creates two new stereocenters.
    # The products are a pair of enantiomers and a meso compound.
    rxn2_products = [
        {"id": "2A", "name": "(2R, 4R)-pentane-2,4-diol", "enantiomer_id": "2B"},
        {"id": "2B", "name": "(2S, 4S)-pentane-2,4-diol", "enantiomer_id": "2A"},
        {"id": "2C", "name": "(2R, 4S)-pentane-2,4-diol (meso)", "enantiomer_id": None}
    ]

    # Combine all products into one list
    all_products = rxn1_products + rxn2_products

    # --- Step 2: Calculate peaks for Chiral HPLC ---
    # Chiral HPLC separates all unique stereoisomers.
    # The number of peaks is simply the total number of unique compounds.
    calculated_chiral_peaks = len(all_products)

    # --- Step 3: Calculate peaks for Normal-Phase HPLC ---
    # Normal-phase (achiral) HPLC separates diastereomers but not enantiomers.
    # We count each diastereomer/meso compound as one peak, and each enantiomeric pair as one peak.
    counted_ids = set()
    calculated_normal_peaks = 0
    for product in all_products:
        if product["id"] in counted_ids:
            continue
        
        # This is a new peak
        calculated_normal_peaks += 1
        
        # Mark this product as counted
        counted_ids.add(product["id"])
        
        # If it has an enantiomer, mark the enantiomer as counted as well,
        # since they will co-elute in an achiral column.
        if product["enantiomer_id"]:
            counted_ids.add(product["enantiomer_id"])

    # --- Step 4: Get the expected values from the provided answer ---
    # The provided answer is <<<C>>>, which corresponds to:
    # "5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC"
    expected_chiral_peaks = 5
    expected_normal_peaks = 4

    # --- Step 5: Verify the answer ---
    if calculated_chiral_peaks == expected_chiral_peaks and calculated_normal_peaks == expected_normal_peaks:
        return "Correct"
    else:
        error_messages = []
        if calculated_chiral_peaks != expected_chiral_peaks:
            error_messages.append(
                f"Incorrect Chiral HPLC peak count. "
                f"Expected: {expected_chiral_peaks}, Calculated: {calculated_chiral_peaks}. "
                f"Reason: There are {len(all_products)} total unique stereoisomers, and a chiral column should resolve all of them."
            )
        if calculated_normal_peaks != expected_normal_peaks:
            error_messages.append(
                f"Incorrect Normal-Phase HPLC peak count. "
                f"Expected: {expected_normal_peaks}, Calculated: {calculated_normal_peaks}. "
                f"Reason: Reaction I produces 2 diastereomers (2 peaks). Reaction II produces an enantiomeric pair (1 peak) and a meso compound (1 peak). Total = 2 + 2 = 4 peaks."
            )
        return "\n".join(error_messages)

# Run the check
result = check_hplc_peaks()
print(result)