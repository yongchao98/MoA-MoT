import re

def check_chemistry_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the organic chemistry question.

    The function verifies two main aspects:
    1. The scientific reasoning is sound and follows the correct logical path.
    2. The final selected option (e.g., <<<A>>>) is consistent with the correct reasoning.
    """

    # Step 1: Define the ground truth based on scientific principles.
    # The jargon "on top of each other" in organic chemistry refers to a failed
    # chromatographic separation (TLC or column chromatography).
    # The principle of this separation method is polarity.
    # Therefore, a failed separation implies the compounds have similar polarities.
    correct_concept = "similar polarities"

    # Step 2: Parse the provided answer to extract its reasoning and final choice.
    reasoning = final_answer_text.lower()
    
    # Extract the final choice letter, e.g., from <<<A>>>
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Incorrect: The final answer does not have the required format '<<<X>>>' at the end."
    
    final_choice_letter = match.group(1)

    # Step 3: Map the chosen letter to its concept based on the answer's own text.
    # The provided answer explicitly lists the options in its "Step 5: Evaluate the Given Options".
    # We can hardcode this mapping as it's defined within the text we are checking.
    options_in_answer = {
        "A": "similar polarities",
        "B": "bonding to each other",
        "C": "similar optical rotations",
        "D": "similar boiling points"
    }
    
    chosen_concept = options_in_answer.get(final_choice_letter)
    if not chosen_concept:
        return f"Incorrect: The final choice '{final_choice_letter}' is not a valid option according to the text's own evaluation section."

    # Step 4: Verify the quality and correctness of the reasoning.
    # A correct explanation must connect the jargon to chromatography and its underlying principle, polarity.
    mentions_chromatography = "chromatography" in reasoning or "tlc" in reasoning
    mentions_polarity = "polarity" in reasoning or "polar" in reasoning
    mentions_separation = "separation" in reasoning or "separate" in reasoning
    
    if not (mentions_chromatography and mentions_polarity and mentions_separation):
        return "Incorrect: The reasoning is flawed. It fails to correctly connect the phrase 'on top of each other' to the concept of chromatographic separation based on polarity."

    # A good explanation should also correctly dismiss the other main distractors.
    dismisses_boiling_point = "distillation" in reasoning and "boiling point" in reasoning
    dismisses_optical_rotation = "optical rotation" in reasoning
    
    if not (dismisses_boiling_point and dismisses_optical_rotation):
        return "Incorrect: The reasoning is incomplete. It does not adequately explain why other options like 'similar boiling points' (distillation) or 'similar optical rotations' are incorrect."

    # Step 5: Perform the final check for consistency and correctness.
    # Does the reasoning correctly identify the right concept?
    if correct_concept not in reasoning:
        return f"Incorrect: The reasoning text does not seem to conclude that '{correct_concept}' is the correct cause."

    # Is the final selected option consistent with the correct concept?
    if chosen_concept == correct_concept:
        return "Correct"
    else:
        return (f"Incorrect: The answer's reasoning correctly points to '{correct_concept}', "
                f"but the final choice is '{final_choice_letter}', which corresponds to '{chosen_concept}'. "
                "There is a contradiction between the reasoning and the final answer.")

# The final answer provided by the user to be checked.
final_answer_to_check = """
Here is a step-by-step analysis of the question to determine the final answer.

**Step 1: Analyze the Context and the Key Phrase**
The conversation takes place between two synthetic organic chemists. This is a crucial detail, as it points to common techniques, problems, and jargon specific to this field. The key phrase, "my compounds are on top of each other," is a colloquialism used to describe a failed experimental outcome. The visual nature of the phrase ("on top of") is a strong clue about the technique being used.

**Step 2: Identify the Most Relevant Laboratory Technique**
In a synthetic organic chemistry lab, the most common technique for both analyzing reaction mixtures and purifying products is **chromatography**, particularly Thin-Layer Chromatography (TLC) and column chromatography.
*   In **TLC**, compounds are separated into distinct spots at different heights on a plate. If two compounds fail to separate, their spots will overlap, literally appearing "on top of each other." This provides a direct and common interpretation of the phrase.
*   In **column chromatography**, compounds move down a column in bands. If the bands do not separate, they co-elute (exit the column together), meaning they were effectively "on top of each other" within the column.
Other techniques like distillation or recrystallization do not fit this visual metaphor as well.

**Step 3: Determine the Underlying Scientific Principle of Chromatography**
The separation of compounds in standard laboratory chromatography (using a polar stationary phase like silica gel) is based almost entirely on the **polarity** of the molecules.
*   The stationary phase (e.g., silica gel) is polar.
*   Polar compounds interact strongly with the stationary phase and move slowly.
*   Less polar compounds are carried along more quickly by the mobile phase (the solvent).
A successful separation relies on the compounds having different polarities, which causes them to move at different rates.

**Step 4: Connect the Problem to the Principle**
A failed chromatographic separation, where compounds are "on top of each other," occurs when the compounds move at the same rate. This happens when their interactions with the stationary and mobile phases are nearly identical. The direct cause for this is that the compounds have very **similar polarities**.

**Step 5: Evaluate the Given Options**
*   **A) The compounds they are working with have similar polarities.** This is the most direct and accurate explanation for a failed chromatographic separation, which is what the chemist's colloquialism describes.
*   **B) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.** This is too general. While these forces are the basis for polarity, the specific problem is that the *net effect* of these interactions (i.e., the overall polarity) is too similar between the different compounds, causing them to interact with the chromatography medium in the same way.
*   **C) The compounds they are working with have similar optical rotations.** Optical rotation is a property of chiral molecules and is not the basis for separation in standard chromatography.
*   **D) The compounds they are working with have similar boiling points.** This would make separation by distillation difficult. However, the phrase "on top of each other" is characteristic of chromatography, not distillation.

**Conclusion**
The chemist is lamenting a failed purification. The phrase "on top of each other" is a classic description of a failed separation by chromatography. The underlying chemical reason for this failure is that the compounds have similar polarities.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(final_answer_to_check)
print(result)