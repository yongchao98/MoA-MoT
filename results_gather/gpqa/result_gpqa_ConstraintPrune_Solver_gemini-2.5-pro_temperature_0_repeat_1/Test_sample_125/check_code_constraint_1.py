import re

def check_correctness_of_hplc_answer(llm_answer_str: str) -> str:
    """
    Checks the correctness of the given answer for the HPLC analysis problem.

    The function models the chemical reactions and HPLC separation principles to
    determine the correct number of peaks and compares it with the provided answer.
    """

    # Step 1: Define the products of each reaction based on stereochemical principles.
    
    # Reaction I: (S)-5-methoxyhexan-3-one + LAH -> 5-methoxyhexan-3-ol
    # Creates a new chiral center at C3, resulting in two diastereomers.
    products_reaction_I = [
        "(3R, 5S)-5-methoxyhexan-3-ol",  # Diastereomer 1
        "(3S, 5S)-5-methoxyhexan-3-ol"   # Diastereomer 2
    ]

    # Reaction II: Pentane-2,4-dione + NaBH4 -> Pentane-2,4-diol
    # Creates two chiral centers (C2, C4), resulting in a pair of enantiomers and a meso compound.
    products_reaction_II = [
        "(2R, 4R)-pentane-2,4-diol",  # Enantiomer A
        "(2S, 4S)-pentane-2,4-diol",  # Enantiomer B (enantiomer of A)
        "(2R, 4S)-pentane-2,4-diol"   # Meso compound (diastereomer of A and B)
    ]

    # The final mixture contains all unique stereoisomers from both reactions.
    final_mixture = products_reaction_I + products_reaction_II

    # Step 2: Calculate the expected number of peaks for each HPLC type.

    # Chiral HPLC separates all unique stereoisomers.
    # The number of peaks is the total number of unique compounds.
    expected_chiral_peaks = len(final_mixture)

    # Normal-phase HPLC separates diastereomers but not enantiomers.
    # Rxn I products: 2 diastereomers -> 2 peaks.
    # Rxn II products: 1 enantiomeric pair (1 peak) + 1 meso compound (1 peak) -> 2 peaks.
    # The products from Rxn I and Rxn II are structurally distinct and will separate.
    # Total peaks = 2 (from Rxn I) + 2 (from Rxn II) = 4.
    expected_normal_phase_peaks = 4

    # Step 3: Parse the LLM's answer and compare it with the calculated correct answer.
    
    match = re.search(r'<<<([A-D])>>>', llm_answer_str)
    if not match:
        return "Failure: Invalid answer format. Expected format like <<<A>>>."
    
    llm_choice = match.group(1)

    # Define the peak counts for each multiple-choice option.
    choices = {
        'A': {'chiral': 3, 'normal': 2},
        'B': {'chiral': 3, 'normal': 3},
        'C': {'chiral': 4, 'normal': 2},
        'D': {'chiral': 5, 'normal': 4}
    }

    llm_peaks = choices[llm_choice]
    
    # Check for correctness and provide detailed feedback if incorrect.
    errors = []
    if llm_peaks['chiral'] != expected_chiral_peaks:
        error_msg = (
            f"Incorrect number of peaks for chiral HPLC. "
            f"The answer suggests {llm_peaks['chiral']} peaks, but the correct number is {expected_chiral_peaks}. "
            f"Reasoning: Chiral HPLC separates all unique stereoisomers. "
            f"Reaction I produces 2 diastereomers. Reaction II produces 3 stereoisomers (an enantiomeric pair and a meso compound). "
            f"This gives a total of 2 + 3 = 5 unique compounds, resulting in 5 peaks."
        )
        errors.append(error_msg)

    if llm_peaks['normal'] != expected_normal_phase_peaks:
        error_msg = (
            f"Incorrect number of peaks for normal-phase HPLC. "
            f"The answer suggests {llm_peaks['normal']} peaks, but the correct number is {expected_normal_phase_peaks}. "
            f"Reasoning: Normal-phase HPLC separates diastereomers but not enantiomers. "
            f"The 2 diastereomers from Reaction I give 2 peaks. "
            f"From Reaction II, the enantiomeric pair co-elutes as 1 peak, and the meso compound gives a separate peak. "
            f"This results in 2 (from Rxn I) + 2 (from Rxn II) = 4 total peaks."
        )
        errors.append(error_msg)

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Given answer from the other LLM
llm_answer = "<<<D>>>"

# Run the check
result = check_correctness_of_hplc_answer(llm_answer)
print(result)