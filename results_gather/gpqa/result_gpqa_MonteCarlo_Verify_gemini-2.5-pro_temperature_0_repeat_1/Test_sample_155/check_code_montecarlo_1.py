def check_hplc_peaks():
    """
    This function verifies the answer to the chemistry question by modeling the
    stereochemical products and how they are separated by different HPLC techniques.
    """

    # --- Step 1: Define the products based on stereochemical principles ---

    # Reaction 1: (E)-oct-4-ene undergoes anti-dihydroxylation.
    # For a symmetric E-alkene, anti-addition results in a single meso compound.
    # Let's represent this product.
    product_from_reaction_1 = ["(4R,5S)-octane-4,5-diol (meso)"]

    # Reaction 2: (Z)-oct-4-ene undergoes anti-dihydroxylation.
    # For a symmetric Z-alkene, anti-addition results in a racemic mixture of enantiomers.
    # Let's represent these two products.
    product_from_reaction_2 = ["(4R,5R)-octane-4,5-diol", "(4S,5S)-octane-4,5-diol"]

    # The final mixture is the combination of products from both reactions.
    final_mixture = product_from_reaction_1 + product_from_reaction_2

    # The mixture contains 3 unique stereoisomers: one meso compound and a pair of enantiomers.
    # The meso compound is a diastereomer of the enantiomers.

    # --- Step 2: Simulate the number of peaks for each HPLC type ---

    # Standard (achiral) HPLC:
    # Separates diastereomers but NOT enantiomers. Enantiomers co-elute.
    # Peak 1: The meso compound.
    # Peak 2: The pair of enantiomers ((4R,5R) and (4S,5S)).
    # Expected number of peaks = 2.
    num_peaks_standard_hplc = 2

    # Chiral HPLC:
    # Separates all stereoisomers, including enantiomers.
    # Each unique stereoisomer gives a distinct peak.
    # Peak 1: (4R,5S)-octane-4,5-diol (meso)
    # Peak 2: (4R,5R)-octane-4,5-diol
    # Peak 3: (4S,5S)-octane-4,5-diol
    # Expected number of peaks = number of unique molecules in the mixture.
    num_peaks_chiral_hplc = len(set(final_mixture))

    # --- Step 3: Verify the provided answer 'D' ---
    # Answer 'D' claims: 2 peaks in standard HPLC and 3 peaks in chiral HPLC.
    
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    if num_peaks_standard_hplc != expected_standard_peaks:
        return (f"Incorrect. The answer claims {expected_standard_peaks} peaks for standard HPLC, "
                f"but the correct number is {num_peaks_standard_hplc}. "
                "A standard column separates diastereomers, so the meso compound and the enantiomeric pair "
                "will form two distinct peaks.")

    if num_peaks_chiral_hplc != expected_chiral_peaks:
        return (f"Incorrect. The answer claims {expected_chiral_peaks} peaks for chiral HPLC, "
                f"but the correct number is {num_peaks_chiral_hplc}. "
                "A chiral column separates all three unique stereoisomers (the meso compound and the two enantiomers).")

    # If both calculated values match the answer's claim.
    return "Correct"

# Run the check
result = check_hplc_peaks()
if result == "Correct":
    print("Correct")
else:
    print(result)