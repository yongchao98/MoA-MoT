import re

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer to a chemistry question.

    The function hardcodes the question, options, and the scientifically correct answer.
    It then parses the provided LLM answer and compares it against the correct one.
    """
    # The original question and its options
    question = "Your overhear two chemists talking to each other as they leave a synthetic organic chemistry lab. One asks the other 'So, how did it go?' The second chemist replies, 'Not well - my compounds are on top of each other.' What is the second chemist most likely referring to?"
    options = {
        "A": "The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.",
        "B": "The compounds they are working with have similar optical rotations.",
        "C": "The compounds they are working with have similar boiling points.",
        "D": "The compounds they are working with have similar polarities."
    }

    # The final answer provided by the LLM to be checked
    llm_answer_text = "<<<D>>>"

    # Scientific Justification for the correct answer
    # 1. Context: The conversation is between synthetic organic chemists, for whom purification is a key task.
    # 2. Jargon: "On top of each other" is a common lab colloquialism for a failed separation. The visual nature of the phrase strongly suggests chromatography (e.g., spots on a TLC plate).
    # 3. Technique: Chromatography is the most common separation technique in this field. It separates compounds based on their differential partitioning between a stationary phase and a mobile phase.
    # 4. Principle: In standard chromatography (e.g., on silica gel), separation is based on polarity.
    # 5. Conclusion: If compounds are "on top of each other," it means they failed to separate via chromatography, which is a direct result of them having very similar polarities.
    correct_answer_key = "D"

    # Parse the LLM's answer
    match = re.search(r'<<<(.+)>>>', llm_answer_text)
    if not match:
        return f"Invalid answer format: The answer '{llm_answer_text}' does not contain '<<<...>>>'."

    llm_answer_key = match.group(1).strip()

    if llm_answer_key not in options:
        return f"Invalid option: The answer '{llm_answer_key}' is not one of the valid options {list(options.keys())}."

    # Check if the LLM's answer is correct
    if llm_answer_key == correct_answer_key:
        return "Correct"
    else:
        # Explain why the LLM's answer is incorrect
        incorrect_reason = ""
        if llm_answer_key == "A":
            incorrect_reason = "This is too general. While these forces are the basis for polarity, the direct cause of the failed separation is the similar *net effect* of these forces, which is best described as 'similar polarities'."
        elif llm_answer_key == "B":
            incorrect_reason = "Similar optical rotation is a property of chiral molecules and is not the basis for separation in standard chromatography."
        elif llm_answer_key == "C":
            incorrect_reason = "Similar boiling points would be a problem for separation by distillation, not chromatography. The phrase 'on top of each other' is characteristic of a failed chromatographic separation."

        return (f"Incorrect. The provided answer is {llm_answer_key}, but the correct answer is {correct_answer_key}.\n"
                f"Reason: The phrase 'on top of each other' is a common lab term for a failed chromatographic separation. This type of separation relies on differences in polarity. Therefore, a failed separation means the compounds have similar polarities. "
                f"The chosen option ({llm_answer_key}) is wrong because: {incorrect_reason}")

# Execute the check and print the result
result = check_chemistry_answer()
print(result)