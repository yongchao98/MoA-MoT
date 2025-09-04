def check_hplc_stereoisomer_separation():
    """
    This function verifies the number of HPLC peaks for a mixture of diol stereoisomers.

    It models the following chemical principles:
    1.  Reaction 1: (E)-oct-4-ene + anti-dihydroxylation -> meso-octane-4,5-diol.
    2.  Reaction 2: (Z)-oct-4-ene + anti-dihydroxylation -> racemic mixture of (4R,5R)- and (4S,5S)-octane-4,5-diol.
    3.  Standard (achiral) HPLC separates diastereomers but not enantiomers.
    4.  Chiral HPLC separates all stereoisomers (diastereomers and enantiomers).
    """

    # Step 1: Define the stereoisomers in the final mixture.
    # We represent each unique molecule with a unique name.
    # The meso compound is a single, achiral stereoisomer.
    meso_diol = "meso-octane-4,5-diol"
    
    # The racemic mixture contains two distinct molecules that are enantiomers of each other.
    enantiomer_A = "(4R,5R)-octane-4,5-diol"
    enantiomer_B = "(4S,5S)-octane-4,5-diol"

    # The final mixture contains all three stereoisomers.
    product_mixture = [meso_diol, enantiomer_A, enantiomer_B]

    # Step 2: Simulate Standard (Achiral) HPLC.
    # On an achiral column, enantiomers co-elute (appear as one peak).
    # Diastereomers elute separately.
    # The meso-diol is a diastereomer to the enantiomeric pair.
    # We can model this by creating "elution groups".
    
    # Group 1: The meso compound.
    # Group 2: The enantiomeric pair, which co-elutes.
    achiral_elution_groups = {
        "meso_group": frozenset([meso_diol]),
        "enantiomeric_pair_group": frozenset([enantiomer_A, enantiomer_B])
    }
    
    calculated_standard_peaks = len(achiral_elution_groups)

    # Step 3: Simulate Chiral HPLC.
    # On a chiral column, all unique stereoisomers are separated.
    # The number of peaks is the number of unique molecules in the mixture.
    calculated_chiral_peaks = len(set(product_mixture))

    # Step 4: Check against the provided answer (B: 2 peaks standard, 3 peaks chiral).
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    if (calculated_standard_peaks == expected_standard_peaks and
            calculated_chiral_peaks == expected_chiral_peaks):
        return "Correct"
    else:
        error_messages = []
        if calculated_standard_peaks != expected_standard_peaks:
            error_messages.append(
                f"Incorrect peak count for standard HPLC. "
                f"Expected: {expected_standard_peaks}, but logic dictates {calculated_standard_peaks}. "
                f"The mixture contains one meso compound and one racemic pair. These are diastereomers, resulting in 2 peaks on an achiral column."
            )
        if calculated_chiral_peaks != expected_chiral_peaks:
            error_messages.append(
                f"Incorrect peak count for chiral HPLC. "
                f"Expected: {expected_chiral_peaks}, but logic dictates {calculated_chiral_peaks}. "
                f"The mixture contains three distinct stereoisomers (meso, R,R-enantiomer, S,S-enantiomer), which are all resolved on a chiral column, resulting in 3 peaks."
            )
        return "\n".join(error_messages)

# Run the check
result = check_hplc_stereoisomer_separation()
print(result)