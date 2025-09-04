def check_answer():
    """
    This function checks the correctness of the provided answer by simulating the
    chemical reactions and chromatographic separations based on stereochemical rules.
    """

    # Step 1 & 2: Define stereochemical rules and determine reaction products
    # We represent the products as unique strings for simplicity.
    # A racemic mixture is represented as a list of its two enantiomers.
    # A meso compound is a single entity.
    
    products = {
        "reaction_1": [],
        "reaction_2": []
    }

    # Reaction 1: (E)-oct-4-ene (trans) undergoes anti-dihydroxylation
    # Rule: trans-alkene + anti-addition -> meso compound
    products["reaction_1"].append("meso-octane-4,5-diol")

    # Reaction 2: (Z)-oct-4-ene (cis) undergoes anti-dihydroxylation
    # Rule: cis-alkene + anti-addition -> racemic mixture
    products["reaction_2"].append("(4R,5R)-octane-4,5-diol")
    products["reaction_2"].append("(4S,5S)-octane-4,5-diol")

    # Step 3: Combine the products into a final mixture
    final_mixture = products["reaction_1"] + products["reaction_2"]
    
    # The final mixture contains 3 distinct stereoisomers:
    # ['meso-octane-4,5-diol', '(4R,5R)-octane-4,5-diol', '(4S,5S)-octane-4,5-diol']
    
    # Step 4: Predict the chromatogram on a standard (achiral) HPLC column
    # On an achiral column, enantiomers are indistinguishable. We can group them.
    # A meso compound is a diastereomer to the enantiomeric pair.
    achiral_groups = set()
    for compound in final_mixture:
        if "meso" in compound:
            achiral_groups.add("meso_group")
        else:
            # Both (R,R) and (S,S) enantiomers belong to the same achiral group
            achiral_groups.add("enantiomeric_pair_group")
            
    standard_hplc_peaks = len(achiral_groups)

    # Step 5: Predict the chromatogram on a chiral HPLC column
    # A chiral column separates all unique stereoisomers.
    # We can simply count the number of unique strings in our mixture list.
    chiral_hplc_peaks = len(set(final_mixture))

    # The question's correct option is D, which corresponds to:
    # 2 peaks in standard HPLC and 3 peaks in chiral HPLC
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    # Check if the calculated results match the expected results from the correct answer.
    if standard_hplc_peaks != expected_standard_peaks:
        return (f"Incorrect. The number of peaks in standard (achiral) HPLC is wrong. "
                f"Expected {expected_standard_peaks} peaks, but the analysis leads to {standard_hplc_peaks} peaks. "
                f"An achiral column should separate the meso compound from the enantiomeric pair, "
                f"but the two enantiomers will co-elute, resulting in 2 peaks total.")

    if chiral_hplc_peaks != expected_chiral_peaks:
        return (f"Incorrect. The number of peaks in chiral HPLC is wrong. "
                f"Expected {expected_chiral_peaks} peaks, but the analysis leads to {chiral_hplc_peaks} peaks. "
                f"A chiral column should resolve all three distinct stereoisomers (the meso compound and the two enantiomers), "
                f"resulting in 3 peaks total.")

    return "Correct"

# Run the check
result = check_answer()
print(result)