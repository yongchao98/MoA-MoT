import re

def check_chemistry_answer():
    """
    Checks the correctness of the LLM's answer to the organic chemistry question.
    """
    # The question's options as provided in the final answer block.
    options = {
        "A": "The compounds they are working with have similar boiling points.",
        "B": "The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.",
        "C": "The compounds they are working with have similar optical rotations.",
        "D": "The compounds they are working with have similar polarities."
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_text = """
    Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

    **1. Deconstruct the Question and Context**
    *   **Setting:** The conversation takes place between two synthetic organic chemists leaving a lab. This context is critical as it points towards common techniques, problems, and jargon specific to this field.
    *   **Core Problem:** One chemist's experiment went "Not well."
    *   **The Jargon:** The reason given is "my compounds are on top of each other." This is a colloquialism, a piece of lab slang, not a formal scientific term. The phrase is highly visual, implying a failure to achieve physical separation.

    **2. Identify the Most Relevant Laboratory Technique**
    The primary goal after a synthesis reaction is purification—separating the desired product from unreacted starting materials, byproducts, and catalysts. The candidate answers overwhelmingly agree that the phrase "on top of each other" points to a specific technique:
    *   **Chromatography:** This is the workhorse technique for separation and analysis in organic chemistry. In Thin-Layer Chromatography (TLC), compounds are visualized as spots at different heights on a plate. If two compounds fail to separate, their spots will overlap, literally appearing "on top of each other." This provides a perfect visual match for the slang term. In column chromatography, the same failure results in the compounds moving together in the same band and eluting from the column at the same time (co-elution).

    **3. Determine the Underlying Scientific Principle**
    The candidate answers correctly identify that the separation in standard laboratory chromatography (using a polar stationary phase like silica gel) is based almost entirely on the **polarity** of the compounds.
    *   The stationary phase (e.g., silica gel) is polar.
    *   Polar compounds interact strongly with the stationary phase and move slowly.
    *   Less polar compounds are carried along more quickly by the less-polar mobile phase (the solvent).
    *   Successful separation depends on the compounds having different polarities.

    **4. Connect the Problem to the Principle**
    If two or more compounds in a mixture have very **similar polarities**, they will interact with the stationary and mobile phases in a nearly identical way. This causes them to travel at the same rate during chromatography, leading to a failed separation. This is the direct cause of the "on top of each other" problem.

    **5. Evaluate and Eliminate the Other Options**
    *   **A) The compounds they are working with have similar boiling points:** This would be a problem for separation by **distillation**. However, the phrase "on top of each other" is the classic jargon for a failed chromatographic separation, not a failed distillation.
    *   **B) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions:** This statement is too general and less precise. While these forces are the fundamental basis for polarity, the direct operational problem is that the *net effect* of these forces (the overall polarity) is too similar for each compound. The problem isn't that they are stuck to each other, but that they are both sticking to the chromatography medium in the same way.
    *   **C) The compounds they are working with have similar optical rotations:** This is an irrelevant property for separation in standard (achiral) chromatography. Enantiomers, for example, have identical polarities and are inseparable by this method, but their optical rotation is a characteristic, not the cause of the separation failure.

    **Conclusion**
    The phrase "my compounds are on top of each other" is a standard colloquialism in organic chemistry for a failed separation using chromatography. The direct chemical reason for this failure is that the compounds have nearly identical polarities, causing them to behave the same way in the chromatographic system.

    <<<D>>>
    """

    # --- Logic for checking the answer ---

    # 1. Establish the correct scientific principle based on the question's context.
    # The phrase "on top of each other" is a classic reference to a failed chromatographic separation.
    # The principle of chromatographic separation is polarity.
    # Therefore, the problem is "similar polarities".
    correct_concept = "similar polarities"

    # 2. Find which option corresponds to the correct concept.
    correct_option_letter = None
    for letter, text in options.items():
        if "similar polarities" in text:
            correct_option_letter = letter
            break
    
    if not correct_option_letter:
        return "Error in checker: The correct concept 'similar polarities' was not found in the options."

    # 3. Extract the letter from the LLM's final answer.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The answer format is incorrect. It should end with '<<<X>>>'."
    
    llm_selected_option = match.group(1)

    # 4. Compare the LLM's answer with the correct answer.
    if llm_selected_option == correct_option_letter:
        # Additionally, check if the reasoning is sound.
        reasoning = llm_answer_text.lower()
        if "chromatography" in reasoning and "polarit" in reasoning: # "polarit" for polarity/polarities
            return "Correct"
        else:
            return f"The selected option <<<{llm_selected_option}>>> is correct, but the reasoning is flawed or incomplete. It fails to mention both chromatography and polarity."
    else:
        return (f"Incorrect. The selected answer is <<<{llm_selected_option}>>>, but the correct answer is <<<{correct_option_letter}>>>. "
                f"The phrase 'on top of each other' refers to a failed chromatographic separation, which occurs when compounds have similar polarities (Option {correct_option_letter}). "
                f"The selected option refers to '{options[llm_selected_option]}', which is incorrect in this context.")

# Run the check and print the result.
result = check_chemistry_answer()
print(result)