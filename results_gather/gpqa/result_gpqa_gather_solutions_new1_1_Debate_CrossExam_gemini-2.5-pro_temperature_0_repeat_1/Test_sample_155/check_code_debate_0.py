def check_answer():
    """
    This function checks the correctness of the answer by modeling the chemical principles involved.
    It simulates the products of the two reactions and then determines the number of peaks
    expected in standard and chiral HPLC based on the stereochemical relationships of the products.
    """

    # Step 1: Define the stereochemical rules for the reaction.
    # The reaction is mCPBA followed by aqueous acid, which is an anti-dihydroxylation.
    # Rule 1: Anti-dihydroxylation of a trans-alkene ((E)-alkene) yields a meso compound.
    # Rule 2: Anti-dihydroxylation of a cis-alkene ((Z)-alkene) yields a racemic mixture of enantiomers.

    # Step 2: Simulate the products of each reaction.
    # Reaction 1 with (E)-oct-4-ene produces one meso compound.
    products_reaction_1 = {'meso-octane-4,5-diol'}

    # Reaction 2 with (Z)-oct-4-ene produces a pair of enantiomers.
    products_reaction_2 = {'(4R,5R)-octane-4,5-diol', '(4S,5S)-octane-4,5-diol'}
    
    # The chemist combines the products.
    final_mixture = products_reaction_1.union(products_reaction_2)

    # Step 3: Simulate Standard (Achiral) HPLC.
    # Principle: Standard HPLC separates diastereomers but not enantiomers.
    # The mixture contains one meso compound and one pair of enantiomers.
    # The meso compound is a diastereomer of the enantiomeric pair.
    # Therefore, the meso compound will form one peak, and the two enantiomers will co-elute as a second peak.
    
    # Count the number of diastereomeric groups.
    # Group 1: The meso compound.
    # Group 2: The enantiomeric pair.
    num_standard_peaks = 2

    # Step 4: Simulate Chiral HPLC.
    # Principle: Chiral HPLC separates all unique stereoisomers (both diastereomers and enantiomers).
    # The number of peaks is equal to the number of unique stereoisomers in the mixture.
    num_chiral_peaks = len(final_mixture)

    # Step 5: Check the calculated results against the expected answer.
    # The provided answer states there should be 2 peaks in standard HPLC and 3 peaks in chiral HPLC (Option D).
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    # Verify the number of unique stereoisomers.
    if len(final_mixture) != 3:
        return (f"Incorrect: The analysis of the product mixture is flawed. "
                f"There should be 3 unique stereoisomers (1 meso, 2 enantiomers), but the logic implies {len(final_mixture)}.")

    # Verify the standard HPLC peak count.
    if num_standard_peaks != expected_standard_peaks:
        return (f"Incorrect: The prediction for standard HPLC is wrong. "
                f"Expected {expected_standard_peaks} peaks, but the analysis leads to {num_standard_peaks} peaks. "
                f"A standard column separates diastereomers, so the meso compound (peak 1) should be separated from the co-eluting enantiomeric pair (peak 2).")

    # Verify the chiral HPLC peak count.
    if num_chiral_peaks != expected_chiral_peaks:
        return (f"Incorrect: The prediction for chiral HPLC is wrong. "
                f"Expected {expected_chiral_peaks} peaks, but the analysis leads to {num_chiral_peaks} peaks. "
                f"A chiral column should resolve all three unique stereoisomers.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)