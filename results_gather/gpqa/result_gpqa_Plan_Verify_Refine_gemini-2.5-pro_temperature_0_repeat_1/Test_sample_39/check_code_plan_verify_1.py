import re

def check_answer(llm_answer_text):
    """
    Checks the correctness of the answer to the chemistry question.

    The function evaluates the provided answer based on common terminology and principles
    in synthetic organic chemistry.

    Args:
        llm_answer_text (str): The full text of the LLM's response, including the final answer.

    Returns:
        str: "Correct" if the answer is correct, or a string explaining the error.
    """
    # Extract the final answer choice (e.g., 'A', 'B', 'C', 'D')
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Invalid answer format. The answer should be in the format <<<X>>> where X is A, B, C, or D."

    llm_choice = match.group(1)
    correct_choice = 'C'

    if llm_choice == correct_choice:
        return "Correct"
    else:
        # Provide a detailed explanation based on chemical principles
        reasoning = {
            'A': "This is incorrect. Optical rotation is a property measured by a polarimeter and does not describe the physical position of compounds during a separation technique like chromatography. The phrase 'on top of each other' refers to a failure to physically separate compounds.",
            'B': "This is incorrect. While similar boiling points make separation by distillation difficult, the phrase 'on top of each other' is a very common idiom specifically for failed *chromatography* (like TLC or column chromatography), where compounds appear as overlapping spots or elute together. The underlying principle for chromatographic separation is polarity, not boiling point.",
            'C': "This is the correct answer. The phrase 'on top of each other' is a common idiom in organic chemistry labs referring to compounds that fail to separate during chromatography (e.g., on a TLC plate or in a column). This failure to separate is a direct result of the compounds having very similar polarities.",
            'D': "This is incorrect because it's too general and less precise than option C. A compound's polarity is the net result of its non-covalent/van der Waals interactions. The problem isn't the existence of these interactions (they are necessary for chromatography to work at all), but that the *overall polarity* is too similar between the compounds, preventing separation. Therefore, C is the more direct and accurate cause of the problem."
        }
        
        error_message = f"The provided answer '{llm_choice}' is incorrect. "
        error_message += f"The correct answer is '{correct_choice}'.\n"
        error_message += f"Reasoning: {reasoning[correct_choice]}\n"
        error_message += f"Why the chosen option '{llm_choice}' is wrong: {reasoning[llm_choice]}"
        
        return error_message

# The provided response from the other LLM
llm_response = """
My initial reasoning led me to conclude that the phrase "on top of each other" refers to a failed separation via chromatography, which is caused by the compounds having similar polarities. The provided feedback confirms that this conclusion is correct. The new information describes liquid-liquid extraction, where immiscible liquids form layers that are literally on top of each other. However, this is typically the *goal* of an extraction, not a problem. A chemist lamenting that things "did not go well" because their compounds are "on top of each other" is almost certainly referring to the frustrating experience of compounds failing to separate on a TLC plate or during column chromatography because their polarities are too similar.

Therefore, the most likely meaning is that the compounds have similar polarities, preventing their separation.

<<<C>>>
"""

# Run the check
result = check_answer(llm_response)
print(result)