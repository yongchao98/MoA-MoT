import re

def check_answer(question, candidate_answers_text, final_answer_text):
    """
    Checks the correctness of the final answer for the given chemistry question.

    Args:
        question (str): The question text.
        candidate_answers_text (str): A string containing all the candidate answers.
        final_answer_text (str): The text containing the final proposed answer.

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """

    # Extract the final answer from the text, e.g., <<<B>>> -> 'B'
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc."
    
    proposed_answer = match.group(1)

    # Define the correct answer and the reasoning based on established chemical principles.
    # This acts as our "ground truth" for the check.
    correct_answer = 'B'
    
    reasoning = {
        'A': "Incorrect. The phrase 'on top of each other' is a visual metaphor strongly associated with chromatography (e.g., spots on a TLC plate). Similar boiling points would be a problem for separation by distillation, for which different jargon (e.g., 'co-distilling') is typically used.",
        'B': "Correct. The phrase 'on top of each other' is common lab jargon for a failed chromatographic separation. Chromatography is a primary purification technique in organic chemistry that separates compounds based on polarity. A failure to separate, where compounds co-elute or appear as a single spot on a TLC plate, is a direct result of them having very similar polarities.",
        'C': "Incorrect. Similar optical rotations are a property of chiral compounds and are not the basis for separation in standard, routine lab chromatography. This property does not cause compounds to be 'on top of each other' in a chromatographic sense.",
        'D': "Incorrect. While non-covalent/van der Waals interactions are the fundamental forces that give rise to polarity, this answer is too general. The direct, practical reason for the failed separation described is that the *net effect* of these interactions results in the compounds having 'similar polarities'. Option B is more precise and directly addresses the principle of chromatographic separation."
    }

    # Check if the proposed answer matches the correct one.
    if proposed_answer == correct_answer:
        return "Correct"
    else:
        # If incorrect, provide the specific reason why that option is wrong.
        return f"Incorrect. The proposed answer is {proposed_answer}, but the correct answer is {correct_answer}. Reason: {reasoning.get(proposed_answer, 'Invalid option selected.')}"

# --- Provided data from the prompt ---

question = """
Your overhear two chemists talking to each other as they leave a synthetic organic chemistry lab. One asks the other "So, how did it go?" The second chemist replies, "Not well - my compounds are on top of each other."

What is the second chemist most likely referring to?

A) The compounds they are working with have similar boiling points.
B) The compounds they are working with have similar polarities.
C) The compounds they are working with have similar optical rotations.
D) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.
"""

# The candidate answers are not needed for the final check, but are included for context.
candidate_answers_text = """
... (all candidate answers from the prompt) ...
"""

# The final answer provided by the LLM to be checked.
final_answer_text = """
Here is a step-by-step analysis of the provided candidate answers to determine the final, correct answer.

**Step 1: Identify the Core Question and Context**
The question asks for the most likely meaning of the phrase "my compounds are on top of each other," as spoken by a synthetic organic chemist. The context is crucial: this is common laboratory jargon.

**Step 2: Analyze the Reasoning Across All Candidate Answers**
Despite the variety in the final lettered answers (`A`, `B`, `C`, `D`), there is a remarkable and near-unanimous consensus in the *reasoning* provided by the agents.

*   **Consensus Point 1: The Technique is Chromatography.** Almost every agent that provides detailed reasoning (e.g., Answers 1, 2, 3, 4, 5, 7, 8, 9, 12, 13, 15) correctly identifies that the phrase "on top of each other" is a visual metaphor strongly associated with a failed separation via chromatography, specifically Thin-Layer Chromatography (TLC) or column chromatography.
*   **Consensus Point 2: The Visual Metaphor.** The agents explain that on a TLC plate, separated compounds appear as distinct spots at different heights. If the compounds fail to separate, their spots overlap, literally appearing "on top of each other." This makes the phrase a direct description of the experimental result.
*   **Consensus Point 3: The Underlying Principle is Polarity.** The agents correctly state that standard chromatography separates compounds based on their differential polarity. Compounds with similar polarities interact with the stationary phase (e.g., silica gel) and mobile phase (solvent) in a nearly identical way.
*   **Consensus Point 4: The Cause of Failure.** The direct consequence of similar polarities is that the compounds travel at the same rate during chromatography, leading to co-elution (failing to separate). This is the problem the chemist is lamenting.

**Step 3: Evaluate the Multiple-Choice Options Based on the Consensus Reasoning**

*   **A) The compounds they are working with have similar boiling points.** The consensus is that this is incorrect. This would be a problem for separation by *distillation*, and the jargon used ("on top of each other") is not typically used for that technique.
*   **B) The compounds they are working with have similar polarities.** This is the correct answer, as established by the overwhelming consensus. It is the direct cause of the failed chromatographic separation that the phrase describes.
*   **C) The compounds they are working with have similar optical rotations.** The consensus is that this is incorrect. Optical rotation is a property of chiral molecules and is not the basis for separation in standard, routine lab chromatography.
*   **D) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.** The consensus is that this is a less precise answer. While these forces are the fundamental basis of polarity, the direct *problem* causing the failed separation is the *similarity* of the net effect of these forces, which is best described as "similar polarities."

**Step 4: Reconcile the Inconsistent Final Answers**
The candidate answers provide a mix of final letters (`A`, `B`, `C`, `D`). However, the detailed explanations almost universally point to "similar polarities" as the correct concept. This strongly suggests that the options were lettered differently in the prompts given to the various agents. For example, Answer 3 provides excellent reasoning for "similar polarities" but concludes with `<<<A>>>`. Answer 15 does the same but concludes with `<<<D>>>`. The only thing that matters is the quality and consistency of the scientific reasoning. The reasoning overwhelmingly supports the concept of similar polarities. In the original question provided, this concept corresponds to option B.

**Conclusion**
The phrase "on top of each other" is a classic piece of organic chemistry slang for a failed separation by chromatography. This failure occurs when the compounds being separated have very similar polarities.

<<<B>>>
"""

# Run the check
result = check_answer(question, candidate_answers_text, final_answer_text)
print(result)