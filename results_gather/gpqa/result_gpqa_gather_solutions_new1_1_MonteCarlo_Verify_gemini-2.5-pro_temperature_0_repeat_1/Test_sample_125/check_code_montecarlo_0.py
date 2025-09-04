def check_hplc_peaks():
    """
    This function programmatically verifies the number of HPLC peaks
    by modeling the products of the two chemical reactions and applying
    the principles of chiral and normal-phase chromatography.
    """

    # --- Step 1: Model the products of each reaction ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one.
    # The starting material is chiral. A new stereocenter is created.
    # The products are a pair of diastereomers.
    # We represent them as two distinct molecules belonging to constitutional group 'A'.
    products_rxn1 = [
        {'id': 'Rxn1_Diastereomer_1', 'constitutional_group': 'A', 'enantiomeric_pair_id': None},
        {'id': 'Rxn1_Diastereomer_2', 'constitutional_group': 'A', 'enantiomeric_pair_id': None}
    ]
    num_products_rxn1 = len(products_rxn1)

    # Reaction II: Reduction of pentane-2,4-dione.
    # The starting material is achiral. Two new stereocenters are created.
    # The products are a pair of enantiomers and a meso compound.
    # We represent them as three distinct molecules belonging to constitutional group 'B'.
    products_rxn2 = [
        {'id': 'Rxn2_Enantiomer_R', 'constitutional_group': 'B', 'enantiomeric_pair_id': 'pair_B1'},
        {'id': 'Rxn2_Enantiomer_S', 'constitutional_group': 'B', 'enantiomeric_pair_id': 'pair_B1'},
        {'id': 'Rxn2_Meso_Compound', 'constitutional_group': 'B', 'enantiomeric_pair_id': None}
    ]
    num_products_rxn2 = len(products_rxn2)

    # The final mixture contains all products from both reactions.
    combined_products = products_rxn1 + products_rxn2

    # --- Step 2: Calculate peaks for Chiral HPLC ---

    # Chiral HPLC separates all unique stereoisomers (diastereomers and enantiomers).
    # The number of peaks is the total number of unique molecules.
    calculated_chiral_peaks = len(combined_products)

    # --- Step 3: Calculate peaks for Normal-Phase HPLC ---

    # Normal-phase (achiral) HPLC separates constitutional isomers and diastereomers,
    # but does NOT separate enantiomers. Enantiomers elute together as one peak.
    # We count peaks by identifying unique "elution groups".
    normal_phase_elution_groups = set()
    for product in combined_products:
        if product['enantiomeric_pair_id'] is not None:
            # Both enantiomers of a pair map to the same elution group identifier.
            group_id = (product['constitutional_group'], product['enantiomeric_pair_id'])
        else:
            # Each diastereomer or meso compound is its own elution group.
            # We use its unique 'id' to distinguish it from others.
            group_id = (product['constitutional_group'], product['id'])
        normal_phase_elution_groups.add(group_id)
        
    calculated_normal_peaks = len(normal_phase_elution_groups)

    # --- Step 4: Verify the candidate answer ---
    
    # The question's correct answer is B: "5 Peaks in chiral HPLC and 4 peaks in normal-phase HPLC"
    expected_chiral_peaks = 5
    expected_normal_peaks = 4

    # Check if our calculated values match the expected values.
    if calculated_chiral_peaks != expected_chiral_peaks:
        return (f"Incorrect. The number of chiral HPLC peaks is wrong. "
                f"Expected {expected_chiral_peaks} peaks, but the analysis leads to {calculated_chiral_peaks} peaks. "
                f"There are {num_products_rxn1} stereoisomers from Reaction I and {num_products_rxn2} from Reaction II, "
                f"totaling {num_products_rxn1 + num_products_rxn2} unique stereoisomers, which should all be separated by chiral HPLC.")

    if calculated_normal_peaks != expected_normal_peaks:
        return (f"Incorrect. The number of normal-phase HPLC peaks is wrong. "
                f"Expected {expected_normal_peaks} peaks, but the analysis leads to {calculated_normal_peaks} peaks. "
                f"Reasoning: Reaction I produces 2 diastereomers (2 peaks). "
                f"Reaction II produces an enantiomeric pair (1 peak) and a meso compound (1 peak), totaling 2 peaks. "
                f"The total should be 2 + 2 = 4 peaks.")

    # If all checks pass, the logic behind the answer is correct.
    return "Correct"

# Run the check
result = check_hplc_peaks()
print(result)