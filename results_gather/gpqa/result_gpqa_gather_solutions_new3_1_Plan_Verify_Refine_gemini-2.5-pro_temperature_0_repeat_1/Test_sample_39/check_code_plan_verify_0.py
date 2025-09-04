import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer to the chemistry question.

    The question asks for the most likely meaning of the phrase "my compounds are on top of each other"
    in a synthetic organic chemistry lab context.

    The correct answer is A, as this phrase is common jargon for a failed chromatographic separation
    (like TLC or column chromatography) due to compounds having similar polarities.

    Args:
        llm_answer_text: The string containing the LLM's full response, including the final answer.

    Returns:
        "Correct" if the answer is correct, or a string explaining the error.
    """
    # Define the correct answer and the reasoning for each option
    correct_option = 'A'
    reasoning = {
        'A': "Correct. The phrase 'on top of each other' is a common colloquialism in organic chemistry labs. It visually describes the result of a failed chromatographic separation (like on a TLC plate) where compounds with similar polarities fail to separate and their spots overlap.",
        'B': "Incorrect. While non-covalent interactions are the fundamental forces at play in chromatography, this answer is too general. The direct, practical reason for the failed separation is that the *net effect* of these interactions, which is described as polarity, is too similar between the compounds.",
        'C': "Incorrect. Similar boiling points would be a problem for separation by distillation. The jargon for a failed distillation is different (e.g., 'they co-distill'). The phrase 'on top of each other' is a visual metaphor strongly associated with chromatography.",
        'D': "Incorrect. Optical rotation is a property of chiral compounds and is not the basis for separation in standard, common laboratory chromatography. Compounds with similar or different optical rotations do not necessarily co-elute for that reason."
    }

    # Extract the final answer from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    
    if not match:
        return "Constraint not satisfied: The final answer was not provided in the required '<<<X>>>' format."

    proposed_option = match.group(1)

    # Check if the proposed answer is the correct one
    if proposed_option == correct_option:
        return "Correct"
    else:
        # Provide the reason why the selected option is incorrect
        return reasoning.get(proposed_option, f"Invalid option '{proposed_option}' selected.")

# The user-provided response to be checked
llm_response = """
Here is a step-by-step analysis of the provided candidate answers to determine the final answer.

**1. Initial Plan: Analyze the Consensus**
My first step is to read through all 15 candidate answers and identify the common thread or consensus in their reasoning. I will tally which of the final options (A, B, C, or D) each candidate supports and why.

**2. Execution: Tallying and Analyzing the Candidate Answers**
*   **Core Context:** Nearly every candidate answer (e.g., 1, 2, 3, 4, 5, 7, 8, 9, 12, 13, 15) correctly identifies the context as a synthetic organic chemistry lab and the phrase "on top of each other" as a common piece of lab jargon.
*   **Key Technique:** There is an overwhelming consensus that this phrase refers to a failed separation attempt using **chromatography**, specifically Thin-Layer Chromatography (TLC) or column chromatography. The phrase is a direct visual metaphor for what a chemist sees on a TLC plate when compounds fail to separate—their spots are unresolved and overlap, appearing "on top of each other."
*   **Underlying Principle:** The candidates unanimously agree that the reason for this failed chromatographic separation is that the compounds in the mixture have **very similar polarities**. Because polarity is the primary property that governs how compounds move in a standard chromatography system, similar polarities lead to similar movement rates and thus, no separation.
*   **Evaluation of Other Options:** The candidates consistently dismiss the other options for clear reasons:
    *   **Similar boiling points (C):** This is relevant for distillation, not chromatography. The jargon for a failed distillation is different (e.g., "they co-distill").
    *   **Similar optical rotations (D):** This is a property of chiral molecules and is not the basis for separation in common lab chromatography.
    *   **Non-covalent interactions (B):** This is too general. While these forces are the physical basis for polarity, the direct and practical reason for the failed separation is the resulting *similar polarity*, which is a more precise answer.

**3. Refined Plan: Synthesize and Conclude**
The analysis shows a near-unanimous conclusion among the candidate answers. The reasoning is consistent and chemically sound. My final step is to synthesize this information and select the best option from the provided list.

*   **Step 1: Identify the core problem.** The chemist is describing a failed separation.
*   **Step 2: Identify the technique.** The phrase "on top of each other" strongly implies a visual technique like chromatography.
*   **Step 3: Identify the chemical principle.** Chromatographic separation is based on polarity. Therefore, a failed separation means the compounds have similar polarities.
*   **Step 4: Match to the options.** Option (A) directly states this conclusion. The other options are less likely or incorrect for the reasons identified by the candidate answers.

The collective reasoning is robust and points decisively to a single answer.

<<<A>>>
"""

# Run the check
result = check_answer_correctness(llm_response)
print(result)