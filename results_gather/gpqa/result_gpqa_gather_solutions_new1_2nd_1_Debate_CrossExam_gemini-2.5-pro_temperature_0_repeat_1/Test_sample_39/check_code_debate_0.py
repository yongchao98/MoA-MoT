import re

def check_answer_correctness(llm_answer: str) -> str:
    """
    Checks if the provided LLM answer to the chemistry question is correct.

    The function verifies two things:
    1. The final selected option (e.g., <<<C>>>) is correct.
    2. The reasoning provided supports the correct answer and dismisses incorrect ones.
    """
    # Define the ground truth for the question
    correct_option = 'C'
    
    # Key concepts that must be in a correct explanation
    required_concepts = {
        "technique": ["chromatography", "tlc"],
        "principle": ["polarity", "polarities"]
    }
    
    # Concepts that should be correctly dismissed
    dismissed_concepts = {
        "technique": ["distillation"],
        "principle": ["boiling point"]
    }

    # 1. Extract the final choice from the answer
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Error: Could not find the final answer in the required format '<<<X>>>'."
    
    llm_choice = match.group(1)
    
    # 2. Check if the final choice is correct
    if llm_choice != correct_option:
        return (f"Incorrect. The final answer choice was '{llm_choice}', but the correct choice is '{correct_option}'. "
                f"The phrase 'on top of each other' refers to a failed chromatographic separation, which is caused by "
                f"compounds having similar polarities.")

    # 3. Check if the reasoning is sound
    reasoning_text = llm_answer.lower()
    
    # Check for presence of required concepts
    if not any(keyword in reasoning_text for keyword in required_concepts["technique"]):
        return ("Incorrect. Although the final choice is correct, the reasoning is flawed. "
                "It fails to mention the key technique of 'chromatography' or 'TLC'.")
                
    if not any(keyword in reasoning_text for keyword in required_concepts["principle"]):
        return ("Incorrect. Although the final choice is correct, the reasoning is flawed. "
                "It fails to mention the key scientific principle of 'polarity'.")

    # Check if the reasoning correctly dismisses the alternative (distillation)
    if not all(keyword in reasoning_text for keyword in dismissed_concepts["technique"] + dismissed_concepts["principle"]):
        return ("Incorrect. The reasoning is incomplete. A full explanation should also address why "
                "other options, such as 'similar boiling points' (related to 'distillation'), are incorrect.")

    # If all checks pass, the answer is correct.
    return "Correct"

# The final answer provided by the user to be checked
final_answer_from_user = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

### Step 1: Analyze the Context and the Core Jargon
The question is set in a synthetic organic chemistry lab. The key to solving the problem is correctly interpreting the lab jargon "my compounds are on top of each other." All 15 candidate answers correctly identify this as a colloquialism describing a negative outcome in a laboratory procedure. The phrase is inherently visual, suggesting a technique where results are observed spatially.

### Step 2: Identify the Relevant Laboratory Technique
There is a unanimous consensus across all 15 candidate answers that the phrase "on top of each other" points directly to the technique of **chromatography**, most likely Thin-Layer Chromatography (TLC). In TLC, chemical compounds in a mixture are separated on a plate and visualized as spots at different heights. If the separation fails, the spots for different compounds overlap, literally appearing "on top of each other." This is a much better fit for the jargon than other techniques like distillation, where one might say compounds "co-distill."

### Step 3: Determine the Scientific Principle Behind the Technique
All candidate answers agree that the primary principle governing separation in common laboratory chromatography (using a silica gel or alumina stationary phase) is **polarity**. The technique works by exploiting the different affinities of compounds for the polar stationary phase and the less-polar mobile phase (the solvent). Compounds with different polarities travel at different rates, which allows them to be separated.

### Step 4: Connect the Problem to the Scientific Principle
A failed chromatographic separation, as described by the chemist, occurs when different compounds travel at the same rate. The reasoning in all candidate answers converges on the same conclusion: this happens when the compounds have very **similar polarities**. Because their polarities are nearly identical, they interact with the stationary and mobile phases in the same way, leading to no separation.

### Step 5: Evaluate the Multiple-Choice Options Based on the Consensus Analysis
Based on the unanimous reasoning of the candidates, we can evaluate the original options:

*   **A) The compounds they are working with have similar optical rotations.** This is incorrect. Optical rotation is a property of chiral molecules and is not the basis for separation in the most common forms of chromatography.
*   **B) The compounds they are working with are bonding to each other through non-covalent/van der Waals interactions.** This is an incorrect choice because it is too general. While these forces are the fundamental physical basis for polarity, the direct, practical problem the chemist is describing is the *consequence* of these forces being too similar. The macroscopic property that governs the separation is polarity, making "similar polarities" the more precise and correct answer.
*   **C) The compounds they are working with have similar polarities.** This is the correct answer. It is the direct cause of a failed chromatographic separation, which is the scenario perfectly described by the chemist's jargon. The reasoning of all 15 candidate answers supports this conclusion.
*   **D) The compounds they are working with have similar boiling points.** This is incorrect. This property is relevant for separation by distillation, and the jargon "on top of each other" is not characteristic of a failed distillation.

### Final Conclusion
Despite the fact that many candidate answers selected an incorrect final letter (likely due to randomization of options during their generation), their underlying reasoning is remarkably consistent and unanimous. They all correctly identify that the chemist is describing a failed chromatographic separation and that the direct cause of this failure is that the compounds have very similar polarities.

<<<C>>>
"""

# Run the check
result = check_answer_correctness(final_answer_from_user)
print(result)