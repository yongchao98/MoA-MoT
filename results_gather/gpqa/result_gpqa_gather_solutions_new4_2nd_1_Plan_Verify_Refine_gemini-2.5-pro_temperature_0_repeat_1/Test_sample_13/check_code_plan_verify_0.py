import re

def check_correctness(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the ChIP-seq question.

    The question asks where ChIP-seq peaks for the IKAROS transcription factor,
    present with PFA fixation, would disappear when using PFA+DSG fixation.

    Args:
        final_answer_text (str): The final answer string, expected in the format "<<<X>>>".

    Returns:
        str: "Correct" if the answer is right, otherwise a string explaining the error.
    """
    # Extract the letter from the answer format, e.g., 'B' from '<<<B>>>'
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format: '{final_answer_text}'. Expected format is <<<X>>>."
    
    selected_answer = match.group(1)
    
    correct_answer = 'B'
    
    # Define the scientific rationale for each option
    reasoning = {
        'A': "Incorrect. ChIP-seq peaks are non-random signals. Their systematic disappearance points to a specific biochemical cause, not a random process.",
        'B': "This is the correct answer. The disappearance of peaks with the addition of DSG (a protein-protein cross-linker) points to a technical artifact like epitope masking or complex insolubility. These artifacts are most severe in regions with the highest local protein density. For a transcription factor like IKAROS, these regions are the active promoters and enhancers where it assembles large regulatory complexes.",
        'C': "Incorrect. While IKAROS can bind to repeats, the most direct explanation for signal *loss* with a *stronger* cross-linker is a biochemical failure (epitope masking/insolubility) driven by the high protein density found at active regulatory sites, which is the hallmark of promoters and enhancers.",
        'D': "Incorrect. This option is too general. While active enhancers can be located in introns, the key feature causing the artifact is the high protein density of the *active regulatory element* itself (promoter or enhancer), not simply its location within an intron."
    }

    if selected_answer == correct_answer:
        return "Correct"
    else:
        return reasoning.get(selected_answer, f"Invalid option '{selected_answer}' selected.")

# The final answer provided was <<<B>>>.
# We will use this as input for our checking function.
final_answer = "<<<B>>>"
result = check_correctness(final_answer)
# The code would return "Correct" for this input.
# If the input were, for example, "<<<C>>>", it would return the explanation for why C is incorrect.