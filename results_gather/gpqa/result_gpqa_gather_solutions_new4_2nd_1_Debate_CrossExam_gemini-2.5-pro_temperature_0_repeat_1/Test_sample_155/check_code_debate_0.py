def check_chemistry_answer():
    """
    This function checks the correctness of the answer by simulating the chemical reactions
    and chromatographic separations based on established stereochemical and analytical principles.
    """

    # --- Step 1 & 2: Define the reaction outcomes based on stereochemical rules ---

    # Let's represent the distinct stereoisomers with unique identifiers.
    MESO_COMPOUND = "meso-octane-4,5-diol"
    RR_ENANTIOMER = "(4R,5R)-octane-4,5-diol"
    SS_ENANTIOMER = "(4S,5S)-octane-4,5-diol"

    # Rule: Anti-dihydroxylation of a trans-alkene gives a meso compound.
    products_reaction_1 = [MESO_COMPOUND]

    # Rule: Anti-dihydroxylation of a cis-alkene gives a racemic mixture.
    products_reaction_2 = [RR_ENANTIOMER, SS_ENANTIOMER]

    # --- Step 3: Combine the products ---
    final_mixture = products_reaction_1 + products_reaction_2
    
    # The final mixture should contain three distinct stereoisomers.
    if len(set(final_mixture)) != 3:
        return f"Constraint check failed: The combined mixture should contain 3 distinct stereoisomers, but the simulation resulted in {len(set(final_mixture))}."

    # --- Step 4: Simulate the HPLC separations ---

    # Standard (achiral) HPLC simulation
    # Enantiomers co-elute. We can simulate this by treating them as a single entity.
    # The meso compound is a diastereomer to the enantiomeric pair and will be separated.
    standard_hplc_peaks = 0
    if MESO_COMPOUND in final_mixture:
        standard_hplc_peaks += 1  # Peak for the meso compound
    if RR_ENANTIOMER in final_mixture or SS_ENANTIOMER in final_mixture:
        standard_hplc_peaks += 1  # One peak for the co-eluting enantiomeric pair
    
    # Chiral HPLC simulation
    # A chiral column separates all unique stereoisomers.
    chiral_hplc_peaks = len(set(final_mixture))

    # --- Step 5: Check the final answer against the simulation ---
    
    # The provided answer states there are 2 peaks in standard HPLC and 3 in chiral HPLC.
    expected_standard_peaks = 2
    expected_chiral_peaks = 3

    if standard_hplc_peaks != expected_standard_peaks:
        return (f"Incorrect analysis for standard HPLC. "
                f"The analysis should yield {expected_standard_peaks} peaks, but the simulation calculated {standard_hplc_peaks}. "
                f"Reason: The meso compound is one peak, and the co-eluting enantiomeric pair is the second peak.")

    if chiral_hplc_peaks != expected_chiral_peaks:
        return (f"Incorrect analysis for chiral HPLC. "
                f"The analysis should yield {expected_chiral_peaks} peaks, but the simulation calculated {chiral_hplc_peaks}. "
                f"Reason: A chiral column separates all three distinct stereoisomers (meso, R,R, and S,S).")

    # The provided answer choice is 'D'. Let's verify if 'D' matches our correct conclusion.
    options = {
        'A': "2 peaks in both standard and chiral HPLC",
        'B': "4 peaks in both standard and chiral HPLC",
        'C': "3 peaks in standard HPLC and 4 peaks in chiral HPLC",
        'D': "2 peaks in standard HPLC and 3 peaks in chiral HPLC"
    }
    
    llm_answer_choice = 'D'
    correct_description = f"{expected_standard_peaks} peaks in standard HPLC and {expected_chiral_peaks} peaks in chiral HPLC"

    if options.get(llm_answer_choice) != correct_description:
        return (f"The reasoning is correct, but the final letter choice is wrong. "
                f"The answer should be '{correct_description}', which corresponds to option D. "
                f"The provided answer chose '{llm_answer_choice}', which is correct.")

    return "Correct"

# Run the check
result = check_chemistry_answer()
print(result)