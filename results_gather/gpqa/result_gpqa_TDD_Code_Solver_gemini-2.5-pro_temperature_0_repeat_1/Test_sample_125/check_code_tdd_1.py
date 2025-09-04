import re

def check_correctness_of_llm_answer(llm_answer_text):
    """
    This function checks the correctness of a given answer to the organic chemistry HPLC question.
    It does so by first determining the correct answer based on chemical principles and then
    comparing it to the provided answer choice.
    """

    # Step 1: Theoretical analysis of the chemical reactions to determine the products.
    # Reaction I: Reduction of (S)-5-methoxyhexan-3-one.
    # The starting material is chiral (S at C5). The reduction of the ketone at C3 creates a new stereocenter.
    # This results in two products: (3R, 5S)-5-methoxyhexan-3-ol and (3S, 5S)-5-methoxyhexan-3-ol.
    # These two products are diastereomers.
    # Number of unique stereoisomers from Reaction I = 2.
    
    # Reaction II: Reduction of pentane-2,4-dione.
    # The starting material is achiral. The reduction of both ketones creates two new stereocenters (C2 and C4).
    # This results in three stereoisomers:
    # 1. A pair of enantiomers: (2R,4R)-pentane-2,4-diol and (2S,4S)-pentane-2,4-diol.
    # 2. A meso compound: (2R,4S)-pentane-2,4-diol, which is achiral.
    # The meso compound is a diastereomer of the enantiomers.
    # Number of unique stereoisomers from Reaction II = 3.

    # Total unique stereoisomers in the final mixture = 2 (from Rxn I) + 3 (from Rxn II) = 5.
    
    # Step 2: Calculate the expected number of peaks for each HPLC type based on the products.
    
    # Chiral HPLC can separate all unique stereoisomers, including enantiomers.
    # Therefore, it will resolve all 5 distinct molecules.
    correct_chiral_peaks = 5
    
    # Normal-phase (achiral) HPLC separates compounds based on polarity. It can separate
    # constitutional isomers and diastereomers, but it cannot separate enantiomers.
    # The peaks observed will be:
    # Peak 1: (3R, 5S)-5-methoxyhexan-3-ol
    # Peak 2: (3S, 5S)-5-methoxyhexan-3-ol (diastereomer of Peak 1, different polarity)
    # Peak 3: The co-eluting racemic mixture of (2R,4R)- and (2S,4S)-pentane-2,4-diol
    # Peak 4: The meso-pentane-2,4-diol (diastereomer of the other diols, different polarity)
    correct_normal_peaks = 4

    # Step 3: Define the options given in the question.
    options = {
        "A": {"chiral": 5, "normal": 4},
        "B": {"chiral": 4, "normal": 2},
        "C": {"chiral": 3, "normal": 3},
        "D": {"chiral": 3, "normal": 2},
    }

    # Step 4: Parse the LLM's answer to find which option (A, B, C, D) was chosen.
    # We use a regular expression to find a standalone letter A, B, C, or D in the response.
    match = re.search(r'\b([A-D])\b', llm_answer_text)
    
    if not match:
        # This handles cases where the LLM gives a non-answer, like the example provided.
        return (f"The provided response is not a valid answer to the multiple-choice question. "
                f"It does not select any of the options A, B, C, or D.\n\n"
                f"For reference, the correct answer is A, which corresponds to {correct_chiral_peaks} chiral peaks and {correct_normal_peaks} normal-phase peaks.")

    llm_choice = match.group(1)
    
    # Step 5: Compare the LLM's chosen option with the correct answer.
    chosen_option_peaks = options[llm_choice]
    
    if chosen_option_peaks["chiral"] == correct_chiral_peaks and chosen_option_peaks["normal"] == correct_normal_peaks:
        return "Correct"
    else:
        # Construct a detailed error message if the choice is incorrect.
        reasons = []
        if chosen_option_peaks["chiral"] != correct_chiral_peaks:
            reasons.append(f"the number of chiral HPLC peaks is incorrect (Answer '{llm_choice}' states {chosen_option_peaks['chiral']}, but it should be {correct_chiral_peaks})")
        if chosen_option_peaks["normal"] != correct_normal_peaks:
            reasons.append(f"the number of normal-phase HPLC peaks is incorrect (Answer '{llm_choice}' states {chosen_option_peaks['normal']}, but it should be {correct_normal_peaks})")
        
        error_message = (
            f"Incorrect. The answer '{llm_choice}' is wrong because {' and '.join(reasons)}.\n\n"
            "**Detailed Analysis:**\n"
            "1.  **Reaction I Products**: Reduction of (S)-5-methoxyhexan-3-one creates a new stereocenter, resulting in two **diastereomers**.\n"
            "2.  **Reaction II Products**: Reduction of pentane-2,4-dione creates two stereocenters, resulting in a **racemic pair of enantiomers** and a **meso compound**.\n"
            "3.  **Total Unique Molecules**: The final mixture contains 5 distinct stereoisomers.\n"
            "4.  **Chiral HPLC**: Separates all 5 unique stereoisomers, resulting in **5 peaks**.\n"
            "5.  **Normal-Phase HPLC**: Separates diastereomers but not enantiomers. This results in 4 distinct groups that can be separated: (diastereomer 1), (diastereomer 2), (the enantiomeric pair), and (the meso compound). This gives **4 peaks**.\n\n"
            f"The correct option is A (5 chiral peaks, 4 normal-phase peaks)."
        )
        return error_message

# The provided answer from the other LLM to be checked.
llm_answer = "Excellent! It's great to see that the Test-Driven Development approach led to the correct answer and passed all checks. I am ready for your next question."

# To check the answer, you would run the function like this:
# result = check_correctness_of_llm_answer(llm_answer)
# print(result)
#
# The expected output for the given llm_answer is:
# "The provided response is not a valid answer to the multiple-choice question. It does not select any of the options A, B, C, or D.
#
# For reference, the correct answer is A, which corresponds to 5 chiral peaks and 4 normal-phase peaks."