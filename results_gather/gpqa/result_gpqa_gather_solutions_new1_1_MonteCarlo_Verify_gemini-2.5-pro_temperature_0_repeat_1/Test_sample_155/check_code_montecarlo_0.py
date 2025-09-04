def check_answer():
    """
    Checks the correctness of the answer to the stereochemistry HPLC problem.
    """
    # The final answer provided by the LLM to be checked.
    # The user's prompt provides a final answer of <<<A>>>.
    llm_answer_choice = 'A'

    # --- Step 1: Define the chemical principles and outcomes ---

    # Reaction 1: (E)-oct-4-ene (trans) undergoes anti-dihydroxylation.
    # Rule: trans-alkene + anti-addition -> meso compound.
    # This produces one unique, achiral compound.
    products_r1 = {"meso_diol"}

    # Reaction 2: (Z)-oct-4-ene (cis) undergoes anti-dihydroxylation.
    # Rule: cis-alkene + anti-addition -> racemic mixture.
    # This produces a pair of enantiomers.
    products_r2 = {"(R,R)_diol", "(S,S)_diol"}

    # The final mixture contains all products.
    # The set of unique stereoisomers is the union of the products.
    final_mixture = products_r1.union(products_r2)
    # final_mixture is now {'meso_diol', '(R,R)_diol', '(S,S)_diol'}
    
    # --- Step 2: Simulate the HPLC separations based on chromatographic rules ---

    # Standard (achiral) HPLC: Separates diastereomers, but not enantiomers.
    # The meso compound is a diastereomer of the (R,R)/(S,S) pair.
    # The (R,R) and (S,S) compounds are enantiomers and will co-elute.
    # Therefore, we expect 2 peaks: one for the meso, one for the enantiomeric pair.
    calculated_standard_peaks = 2

    # Chiral HPLC: Separates all stereoisomers (diastereomers and enantiomers).
    # Therefore, we expect a peak for each unique stereoisomer in the mixture.
    calculated_chiral_peaks = len(final_mixture) # This will be 3

    # --- Step 3: Compare the calculated results with the chosen answer option ---

    # Define the outcomes for each multiple-choice option.
    options = {
        'A': {'standard': 2, 'chiral': 3},
        'B': {'standard': 3, 'chiral': 4},
        'C': {'standard': 4, 'chiral': 4},
        'D': {'standard': 2, 'chiral': 2}
    }

    if llm_answer_choice not in options:
        return f"Invalid answer choice '{llm_answer_choice}'. The choice must be A, B, C, or D."

    expected_results = options[llm_answer_choice]

    # Check if the calculated results match the results from the chosen option.
    is_standard_correct = (calculated_standard_peaks == expected_results['standard'])
    is_chiral_correct = (calculated_chiral_peaks == expected_results['chiral'])

    if is_standard_correct and is_chiral_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_standard_correct:
            error_messages.append(
                f"Standard HPLC peak count is wrong. The answer expects {expected_results['standard']} peaks, "
                f"but the correct analysis yields {calculated_standard_peaks} peaks. "
                "Reason: An achiral column separates the meso compound from the enantiomeric pair, but the two enantiomers co-elute, resulting in 2 peaks."
            )
        if not is_chiral_correct:
            error_messages.append(
                f"Chiral HPLC peak count is wrong. The answer expects {expected_results['chiral']} peaks, "
                f"but the correct analysis yields {calculated_chiral_peaks} peaks. "
                "Reason: A chiral column resolves all three distinct stereoisomers (the meso compound and the two enantiomers), resulting in 3 peaks."
            )
        return "\n".join(error_messages)

# Run the check and print the result.
print(check_answer())