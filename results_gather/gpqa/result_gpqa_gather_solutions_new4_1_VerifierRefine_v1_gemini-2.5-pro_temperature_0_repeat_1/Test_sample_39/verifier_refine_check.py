import re

def check_chemistry_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the chemistry jargon question.

    The question asks for the meaning of "my compounds are on top of each other"
    in a synthetic organic chemistry context. The options are:
    A) Similar boiling points.
    B) Bonding to each other through non-covalent/van der Waals interactions.
    C) Similar optical rotations.
    D) Similar polarities.

    This function verifies if the selected answer is D, which is the correct choice.
    """
    # The correct answer key based on expert knowledge of chemistry lab terminology.
    correct_answer_key = 'D'

    # Reasoning for why each option is correct or incorrect.
    reasoning = {
        'A': "Incorrect. The selected answer is A. Similar boiling points are a problem for separation by distillation. The phrase 'on top of each other' is a visual metaphor that strongly suggests a failed chromatographic separation (like on a TLC plate), not distillation.",
        'B': "Incorrect. The selected answer is B. While non-covalent interactions are the fundamental forces at play, 'similar polarities' is the direct, practical cause of the failed separation. Polarity is the macroscopic property that results from the sum of these interactions and governs chromatographic separation.",
        'C': "Incorrect. The selected answer is C. Similar optical rotations relate to chirality. The root cause of separation failure in standard chromatography is the similarity in polarity, not optical rotation.",
        'D': "Correct. The phrase 'on top of each other' is a common colloquialism in organic chemistry for a failed separation using chromatography (like TLC or column chromatography). This technique separates compounds based on polarity. If compounds have very similar polarities, they will not separate and will appear to be 'on top of each other' on a TLC plate or co-elute from a column."
    }

    # Use regex to find the answer in the format <<<X>>>
    match = re.search(r'<<<([A-D])>>>', final_answer_text)

    if not match:
        return "Error: Could not find a final answer in the required format '<<<X>>>' in the provided text."

    submitted_answer_key = match.group(1)

    # Compare the submitted answer with the correct answer.
    if submitted_answer_key == correct_answer_key:
        return "Correct"
    else:
        # Return the specific reason why the submitted answer is wrong.
        return reasoning.get(submitted_answer_key, f"Invalid option '{submitted_answer_key}' was provided.")

# The final answer provided by the LLM to be checked.
llm_final_answer = """
Here is a step-by-step analysis of the question and the provided answers.

**Step 1: Deconstruct the Question and Context**
The question asks for the most likely meaning behind a chemist's statement: "my compounds are on top of each other." The context is crucial: it's a conversation between two synthetic organic chemists. This field heavily relies on the synthesis of new molecules and, just as importantly, their purification and analysis.

**Step 2: Analyze the Colloquial Phrase**
The phrase "on top of each other" is a common piece of lab jargon. It implies a lack of separation between two or more substances in a mixture. The phrasing is visual, suggesting a technique where the results are observed spatially.

**Step 3: Identify the Most Relevant Laboratory Technique**
In synthetic organic chemistry, the workhorse technique for both analyzing reaction progress and purifying products is **chromatography**.
*   **Thin-Layer Chromatography (TLC):** Used for quick analysis. A spot of the mixture is placed on a plate, and a solvent moves up the plate, separating the components. The result is a series of spots at different heights.
*   **Column Chromatography:** Used for purification. The mixture is passed through a column packed with an adsorbent material. Different compounds travel through the column at different rates and are collected separately.

In both TLC and column chromatography, if two compounds are not separated, they will co-locate. On a TLC plate, their spots will overlap, literally appearing "on top of each other." In a column, they will elute at the same time (co-elute). This visual and spatial outcome makes chromatography the most likely context for the phrase.

**Step 4: Determine the Underlying Principle of Chromatographic Separation**
The most common forms of laboratory chromatography (normal-phase on silica gel or alumina) separate compounds based on their **polarity**. The stationary phase (e.g., silica gel) is polar. Polar compounds in the mixture will adhere more strongly to the stationary phase and move more slowly. Less polar compounds will be carried along more quickly by the mobile phase (the solvent).

**Step 5: Connect the Principle to the Problem**
If two or more compounds in the mixture have very **similar polarities**, they will interact with the stationary and mobile phases in almost the exact same way. Consequently, they will travel at the same rate during the chromatographic process. This results in a failed separation, which is visually represented as spots "on top of each other" on a TLC plate.

**Step 6: Evaluate the Answer Choices**
*   **A) The compounds they are working with have similar boiling points.** This is a problem for separation by **distillation**. While a valid purification issue, the phrase "on top of each other" is much more characteristic of the visual result of chromatography.
*   **B) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.** While these forces are the fundamental basis for physical properties like polarity, this statement is too general. The specific problem is not that the compounds are sticking to *each other*, but that their individual interactions with the chromatography medium are too similar due to their similar overall polarity.
*   **C) The compounds they are working with have similar optical rotations.** Optical rotation is a physical property of chiral molecules and is not the basis for separation in standard chromatography. Enantiomers have identical polarities and are inseparable by these methods, but the root cause of the separation failure is the identical polarity, not the optical rotation.
*   **D) The compounds they are working with have similar polarities.** This is the most direct and accurate explanation. Similar polarities directly cause the failed chromatographic separation that is so aptly described by the phrase "on top of each other."

**Conclusion**
The overwhelming consensus from the analysis is that the chemist is referring to a failed separation by chromatography, and the underlying chemical reason for this failure is that the compounds have similar polarities.

<<<D>>>
"""

# Run the check
result = check_chemistry_answer(llm_final_answer)
print(result)