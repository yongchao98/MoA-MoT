def check_chemistry_hplc_answer():
    """
    Checks the correctness of the answer to the stereochemistry and HPLC problem.

    The function simulates the chemical reactions and chromatographic separations
    based on established rules of organic chemistry.
    """

    # --- Step 1 & 2: Determine the products of each reaction ---

    # The reaction sequence (mCPBA followed by H3O+) is an anti-dihydroxylation.
    
    # Rule for Reaction 1: Anti-dihydroxylation of a trans-alkene ((E)-oct-4-ene)
    # results in a single meso compound.
    products_reaction_1 = {"meso-octane-4,5-diol"}
    
    # Rule for Reaction 2: Anti-dihydroxylation of a cis-alkene ((Z)-oct-4-ene)
    # results in a racemic mixture of two enantiomers.
    products_reaction_2 = {"(4R,5R)-octane-4,5-diol", "(4S,5S)-octane-4,5-diol"}
    
    # --- Step 3: Analyze the combined product mixture ---
    
    # The final mixture contains all unique products from both reactions.
    final_mixture = products_reaction_1.union(products_reaction_2)
    # The mixture is: {'meso-octane-4,5-diol', '(4R,5R)-octane-4,5-diol', '(4S,5S)-octane-4,5-diol'}
    
    # --- Step 4: Predict the number of peaks in Standard (Achiral) HPLC ---
    
    # An achiral column separates diastereomers but not enantiomers.
    # We can group the mixture into sets of compounds that will co-elute.
    # The meso compound is a diastereomer to the other two.
    # The (4R,5R) and (4S,5S) compounds are enantiomers and will co-elute.
    achiral_hplc_groups = [
        {"meso-octane-4,5-diol"},  # Peak 1
        {"(4R,5R)-octane-4,5-diol", "(4S,5S)-octane-4,5-diol"}  # Peak 2
    ]
    correct_standard_peaks = len(achiral_hplc_groups) # Should be 2

    # --- Step 5: Predict the number of peaks in Chiral HPLC ---

    # A chiral column separates all distinct stereoisomers (both diastereomers and enantiomers).
    # The number of peaks will be equal to the number of unique stereoisomers in the mixture.
    correct_chiral_peaks = len(final_mixture) # Should be 3

    # --- Step 6: Check the provided final answer ---

    # The provided final answer is <<<D>>>.
    provided_answer_option = "D"
    
    # Define the claims of each multiple-choice option.
    options = {
        "A": {"standard": 4, "chiral": 4},
        "B": {"standard": 2, "chiral": 2},
        "C": {"standard": 3, "chiral": 4},
        "D": {"standard": 2, "chiral": 3}
    }
    
    if provided_answer_option not in options:
        return f"Invalid answer option '{provided_answer_option}' provided."
        
    selected_option_values = options[provided_answer_option]
    
    # Compare the calculated correct values with the values from the selected option.
    is_standard_correct = (correct_standard_peaks == selected_option_values["standard"])
    is_chiral_correct = (correct_chiral_peaks == selected_option_values["chiral"])
    
    if is_standard_correct and is_chiral_correct:
        return "Correct"
    else:
        error_message = ""
        if not is_standard_correct:
            error_message += (f"The number of peaks in standard HPLC is incorrect. "
                              f"The answer claims {selected_option_values['standard']} peaks, but the correct number is {correct_standard_peaks}. "
                              f"A standard (achiral) column separates the meso compound from the enantiomeric pair, resulting in 2 peaks.\n")
        if not is_chiral_correct:
            error_message += (f"The number of peaks in chiral HPLC is incorrect. "
                              f"The answer claims {selected_option_values['chiral']} peaks, but the correct number is {correct_chiral_peaks}. "
                              f"A chiral column resolves all three distinct stereoisomers (the meso compound and the two enantiomers), resulting in 3 peaks.")
        return error_message.strip()

# Execute the check
result = check_chemistry_hplc_answer()
print(result)