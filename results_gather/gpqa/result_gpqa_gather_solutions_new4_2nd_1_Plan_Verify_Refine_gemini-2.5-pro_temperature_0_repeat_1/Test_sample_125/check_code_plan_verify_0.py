import re

def check_hplc_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the HPLC problem.

    The function simulates the chemical reactions and HPLC separation principles
    to determine the correct number of peaks and compares it with the LLM's reasoning and final choice.
    """

    # --- Step 1: Correct Chemical Analysis (Ground Truth) ---

    # Reaction I: (S)-5-methoxyhexan-3-one (chiral) is reduced at C3.
    # The original stereocenter at C5 is unaffected. A new stereocenter is created at C3.
    # This results in two products: (3R, 5S)-... and (3S, 5S)-...
    # These products are diastereomers.
    # Number of unique products from Reaction I = 2.
    num_products_rxn1 = 2

    # Reaction II: Pentane-2,4-dione (achiral) is reduced at C2 and C4.
    # Two new stereocenters are created.
    # This results in three stereoisomers:
    # 1. (2R, 4R)-pentane-2,4-diol
    # 2. (2S, 4S)-pentane-2,4-diol (enantiomer of 1)
    # 3. (2R, 4S)-pentane-2,4-diol (a meso compound, diastereomer of 1 and 2)
    # Number of unique products from Reaction II = 3.
    num_products_rxn2 = 3

    # --- Step 2: Correct HPLC Peak Calculation (Ground Truth) ---

    # Normal-Phase HPLC (achiral column):
    # Separates diastereomers and constitutional isomers, but NOT enantiomers.
    # Peaks from Rxn I: The 2 diastereomers will separate -> 2 peaks.
    # Peaks from Rxn II: The enantiomeric pair co-elutes (1 peak), and the meso compound separates (1 peak) -> 2 peaks.
    # Total normal-phase peaks = 2 (from Rxn I) + 2 (from Rxn II) = 4.
    correct_normal_peaks = 2 + 2

    # Chiral HPLC (chiral column):
    # Separates all unique stereoisomers (diastereomers and enantiomers).
    # We simply count the total number of unique molecules.
    # Total chiral peaks = 2 (from Rxn I) + 3 (from Rxn II) = 5.
    correct_chiral_peaks = num_products_rxn1 + num_products_rxn2

    # --- Step 3: Parse the LLM's Answer ---
    
    # Extract the final letter choice, e.g., <<<B>>>
    final_choice_match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not final_choice_match:
        return "Incorrect: The final answer is not in the required format '<<<X>>>'."
    llm_choice = final_choice_match.group(1)

    # Extract the numbers from the LLM's reasoning text.
    # Look for phrases like "5 peaks in chiral" and "4 peaks in normal-phase"
    chiral_peaks_match = re.search(r'(\d+)\s*peaks\s*in\s*(the)?\s*chiral\s*hplc', llm_answer_text, re.IGNORECASE)
    normal_peaks_match = re.search(r'(\d+)\s*peaks\s*in\s*(the)?\s*normal-phase\s*hplc', llm_answer_text, re.IGNORECASE)

    if not chiral_peaks_match or not normal_peaks_match:
        # Try another common phrasing
        chiral_peaks_match = re.search(r'chiral\s*hplc\D*(\d+)\s*peaks', llm_answer_text, re.IGNORECASE)
        normal_peaks_match = re.search(r'normal-phase\s*hplc\D*(\d+)\s*peaks', llm_answer_text, re.IGNORECASE)

    if not chiral_peaks_match or not normal_peaks_match:
        return "Incorrect: The reasoning does not clearly state the number of peaks for both chiral and normal-phase HPLC. For example, it should contain '5 peaks in chiral HPLC' and '4 peaks in normal-phase HPLC'."

    llm_chiral_peaks = int(chiral_peaks_match.group(1))
    llm_normal_peaks = int(normal_peaks_match.group(1))

    # --- Step 4: Define the options from the question ---
    options = {
        "A": {"chiral": 3, "normal": 2},
        "B": {"chiral": 5, "normal": 4},
        "C": {"chiral": 4, "normal": 2},
        "D": {"chiral": 3, "normal": 3}
    }

    # --- Step 5: Compare and Validate ---

    # Check if the LLM's reasoning is correct
    if llm_chiral_peaks != correct_chiral_peaks:
        return f"Incorrect: The reasoning is flawed. It states there are {llm_chiral_peaks} peaks in chiral HPLC, but the correct number is {correct_chiral_peaks}."
    
    if llm_normal_peaks != correct_normal_peaks:
        return f"Incorrect: The reasoning is flawed. It states there are {llm_normal_peaks} peaks in normal-phase HPLC, but the correct number is {correct_normal_peaks}."

    # Check if the LLM's final choice is consistent with its (now verified) correct reasoning
    chosen_option_values = options.get(llm_choice)
    if chosen_option_values["chiral"] != correct_chiral_peaks or chosen_option_values["normal"] != correct_normal_peaks:
        return f"Incorrect: The final choice '{llm_choice}' corresponds to {chosen_option_values['chiral']} chiral peaks and {chosen_option_values['normal']} normal peaks. This contradicts the correct analysis of {correct_chiral_peaks} chiral and {correct_normal_peaks} normal peaks, which is option 'B'."

    # If all checks pass, the answer is correct.
    return "Correct"
