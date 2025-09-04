def check_correctness():
    """
    This function checks the correctness of the final answer to the chemistry question.
    It does this by programmatically determining the correct number of HPLC peaks
    and comparing it to the provided answer.
    """
    
    # --- Step 1: Define the chemical logic based on the question ---

    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one (a single enantiomer).
    # This creates a new stereocenter at C3, resulting in two diastereomers.
    # e.g., (3R, 5S) and (3S, 5S) products.
    # Number of unique stereoisomers from Reaction I = 2.
    num_products_rxn1 = 2

    # Reaction II: Reduction of pentane-2,4-dione (achiral).
    # This creates two new stereocenters (C2, C4) in a symmetric molecule.
    # This results in one enantiomeric pair ((2R,4R) & (2S,4S)) and one meso compound ((2R,4S)).
    # Number of unique stereoisomers from Reaction II = 3.
    num_products_rxn2 = 3

    # --- Step 2: Calculate the correct number of peaks for each HPLC type ---

    # Normal-Phase HPLC (achiral column):
    # Separates diastereomers, but not enantiomers.
    # Peaks from Rxn I: The 2 diastereomers are separable. -> 2 peaks.
    # Peaks from Rxn II: The enantiomeric pair co-elutes (1 peak). The meso compound is a diastereomer to the pair and separates (1 peak). -> 2 peaks.
    # Total normal peaks = 2 (from Rxn I) + 2 (from Rxn II)
    correct_normal_peaks = 4

    # Chiral HPLC (chiral column):
    # Separates all unique stereoisomers (both diastereomers and enantiomers).
    # Total chiral peaks = Total number of unique stereoisomers.
    # Total chiral peaks = 2 (from Rxn I) + 3 (from Rxn II)
    correct_chiral_peaks = 5

    # --- Step 3: Parse the provided answer and check against the correct values ---
    
    # The final answer from the provided text is <<<A>>>.
    final_answer_letter = 'A'

    options = {
        'A': {'chiral': 5, 'normal': 4},
        'B': {'chiral': 4, 'normal': 2},
        'C': {'chiral': 3, 'normal': 3},
        'D': {'chiral': 3, 'normal': 2}
    }

    if final_answer_letter not in options:
        return f"Invalid answer option '{final_answer_letter}' found. Cannot check correctness."

    selected_option_values = options[final_answer_letter]
    
    is_chiral_correct = selected_option_values['chiral'] == correct_chiral_peaks
    is_normal_correct = selected_option_values['normal'] == correct_normal_peaks

    if is_chiral_correct and is_normal_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_chiral_correct:
            error_messages.append(
                f"the number of peaks in chiral HPLC is incorrect. The correct number is {correct_chiral_peaks}, but the answer claims {selected_option_values['chiral']}"
            )
        if not is_normal_correct:
            error_messages.append(
                f"the number of peaks in normal-phase HPLC is incorrect. The correct number is {correct_normal_peaks}, but the answer claims {selected_option_values['normal']}"
            )
        
        return f"Incorrect because {', and '.join(error_messages)}."

# Execute the check and print the result.
result = check_correctness()
print(result)