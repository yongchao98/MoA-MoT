import re

def check_chemistry_jargon_answer(question, options, final_answer_text):
    """
    Checks the correctness of an answer to a question about organic chemistry lab jargon.

    This function encodes the domain-specific knowledge required to answer the question
    and compares the logic-derived correct answer with the provided one.

    Args:
        question (str): The question text.
        options (dict): A dictionary mapping option letters to their text.
        final_answer_text (str): The final answer provided by the LLM, including the <<<>>> format.

    Returns:
        str: "Correct" if the answer is correct, otherwise a string explaining the error.
    """
    
    # Extract the letter from the final answer format, e.g., "<<<D>>>" -> "D"
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return f"Invalid answer format. Expected '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>', but got: {final_answer_text}"
    
    provided_answer_letter = match.group(1)

    # --- Domain Knowledge & Rule Engine ---
    # 1. Analyze the key phrase and context from the question.
    key_phrase = "on top of each other"
    context = "synthetic organic chemistry lab"

    # 2. Establish rules based on domain knowledge.
    # Rule 1: The phrase "on top of each other" is a visual metaphor used to describe a failed separation.
    # Rule 2: The most common lab technique with a direct visual output of separation is chromatography (TLC, column).
    # Rule 3: The primary principle governing separation in standard chromatography is polarity.
    # Rule 4: Therefore, a failed separation where compounds are "on top of each other" is most directly caused by the compounds having similar polarities.

    # 3. Evaluate each option against the rules.
    correct_option = None
    reasoning = {}

    # Option A: Similar boiling points
    if "boiling points" in options['A']:
        reasoning['A'] = "This relates to separation by distillation. The jargon for a failed distillation is different (e.g., 'co-distillation'). The phrase 'on top of each other' is strongly associated with the visual output of chromatography."

    # Option B: Similar optical rotations
    if "optical rotations" in options['B']:
        reasoning['B'] = "This is a physical property measured by a polarimeter and is not the basis for common separation techniques like TLC or column chromatography."

    # Option C: Non-covalent interactions
    if "non-covalent" in options['C']:
        reasoning['C'] = "This is too general. While these interactions are the fundamental basis for polarity and how chromatography works, 'similar polarities' is the direct, practical cause of the failed separation that a chemist would refer to."

    # Option D: Similar polarities
    if "similar polarities" in options['D']:
        reasoning['D'] = "This is the correct answer. A failed chromatographic separation (where spots on a TLC plate overlap) is directly caused by compounds having similar polarities, making them travel at the same rate."
        correct_option = 'D'

    # --- Verification ---
    if provided_answer_letter == correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is '{provided_answer_letter}', but the correct answer is '{correct_option}'.\n"
                f"Reasoning: The phrase 'on top of each other' is a common colloquialism in organic chemistry for a failed chromatographic separation. "
                f"This type of separation relies on differences in polarity. If compounds have similar polarities, they fail to separate, appearing as overlapping spots on a TLC plate. "
                f"The reason option '{provided_answer_letter}' is wrong is: {reasoning.get(provided_answer_letter, 'It does not correctly interpret the lab jargon.')}")


# --- Data from the problem ---
question_text = "Your overhear two chemists talking to each other as they leave a synthetic organic chemistry lab. One asks the other 'So, how did it go?' The second chemist replies, 'Not well - my compounds are on top of each other.' What is the second chemist most likely referring to?"

options_dict = {
    'A': "The compounds they are working with have similar boiling points.",
    'B': "The compounds they are working with have similar optical rotations.",
    'C': "The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.",
    'D': "The compounds they are working with have similar polarities."
}

# The final answer provided by the LLM analysis
llm_final_answer = "<<<D>>>"

# --- Execute the check ---
result = check_chemistry_jargon_answer(question_text, options_dict, llm_final_answer)
print(result)