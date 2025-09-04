import re

def check_correctness_of_nmr_answer(final_answer_text):
    """
    Checks the correctness of the final answer for the given NMR problem.

    The logic is based on deducing the structure (para-haloacetophenone) and
    predicting its 1H NMR spectrum.
    """
    # Extract the letter from the final answer format, e.g., <<<D>>>
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not parse the final answer. Expected format is <<<X>>>."
    
    llm_answer = match.group(1)

    # --- Define the constraints and expected features ---
    # The deduced structure is para-haloacetophenone.
    # Expected features in the 1H NMR spectrum:
    # 1. A singlet for 3 protons (acetyl group -COCH3).
    # 2. Two doublets, each for 2 protons (para-substituted aromatic ring).
    # 3. Aromatic signals should be in the ~6.5-8.5 ppm range.
    # 4. Acetyl signal should be in the ~2.0-2.7 ppm range.

    options_analysis = {
        'A': {
            "is_correct": False,
            "reason": "This option is incorrect because the signal at 9.9 ppm is characteristic of an aldehyde, not the ketone structure deduced. It also lacks the required 3H singlet for an acetyl group."
        },
        'B': {
            "is_correct": False,
            "reason": "This option is incorrect because it lacks any signals in the typical aromatic region (6.5-8.5 ppm), which is a required feature for an aromatic compound."
        },
        'C': {
            "is_correct": False,
            "reason": "This option is incorrect because a di-substituted benzene ring must have 4 aromatic protons, but this spectrum only shows one (at 6.9 ppm)."
        },
        'D': {
            "is_correct": True,
            "reason": "This option is correct as it perfectly matches the predicted spectrum for para-haloacetophenone: two doublets for the 4 aromatic protons and a 3H singlet for the acetyl group."
        }
    }

    correct_option_key = 'D'

    if llm_answer == correct_option_key:
        return "Correct"
    else:
        return f"Incorrect. The provided answer was '{llm_answer}'. {options_analysis[llm_answer]['reason']} The correct answer is '{correct_option_key}'."

# The final answer provided in the prompt is <<<D>>>.
# We can run the check on this answer.
final_answer_from_llm = "<<<D>>>" 
# print(check_correctness_of_nmr_answer(final_answer_from_llm))