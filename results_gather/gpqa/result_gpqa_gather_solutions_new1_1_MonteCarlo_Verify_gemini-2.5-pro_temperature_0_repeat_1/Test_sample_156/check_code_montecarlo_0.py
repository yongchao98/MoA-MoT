def check_answer(selected_answer):
    """
    Checks the correctness of the selected answer for designing a retrovirus diagnostic kit.

    The function evaluates the answer based on three core constraints:
    1. The method must be appropriate for a retrovirus (RNA genome).
    2. The kit must be a "molecular" diagnostic kit (detects nucleic acids).
    3. The method must be suitable for "quick detection" (fast and for early diagnosis).
    """

    # Properties of each option based on scientific principles
    options_analysis = {
        'A': {
            'is_correct': False,
            'reason': "The initial identification method is flawed. One cannot design a specific molecular test based on non-specific patient symptoms."
        },
        'B': {
            'is_correct': False,
            'reason': "The identification method is incorrect. A retrovirus has an RNA genome, so direct DNA sequencing is not the correct first step. It requires conversion to cDNA first."
        },
        'C': {
            'is_correct': True,
            'reason': "This option correctly identifies the need for cDNA sequencing for an RNA virus and proposes using real-time PCR, which is the gold standard for quick and accurate molecular diagnosis."
        },
        'D': {
            'is_correct': False,
            'reason': "This method is incorrect for two reasons. First, an ELISA kit that detects antibodies is an immunological test, not a 'molecular' one. Second, IgG antibodies are markers for a later-stage or past infection, making them unsuitable for 'quick detection' in an early outbreak."
        }
    }

    if selected_answer not in options_analysis:
        return f"Invalid option '{selected_answer}'. Please choose from A, B, C, D."

    result = options_analysis[selected_answer]

    if result['is_correct']:
        # Verify that the provided answer matches the correct option
        if selected_answer == 'C':
            return "Correct"
        else:
            # This case handles if the logic determined another option was correct, which would be an error in the logic itself.
            return f"Logic error: The checker determined {selected_answer} is correct, but the established correct answer is C."
    else:
        # The selected answer is incorrect, return the reason.
        return f"Incorrect. The answer '{selected_answer}' is wrong because: {result['reason']}"

# The final answer provided by the LLM is 'C'.
# We will run the check on this answer.
final_answer_to_check = 'C'
print(check_answer(final_answer_to_check))