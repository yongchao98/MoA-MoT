import re

def check_correctness(question, llm_answer):
    """
    Checks the correctness of the LLM's answer for the given organic chemistry question.

    The function verifies two main chemical principles:
    1. The identity of the reactant (A) for an enamine synthesis.
    2. The suitability of the catalyst (B) for an acid-catalyzed dehydration.
    """

    # --- Step 1: Define the chemical principles and correct components ---

    # Principle 1: The reaction is an enamine synthesis from a secondary amine (3-methylpyrrolidine)
    # and a carbonyl compound. The product structure, 1-(cyclohexylidenemethyl)-3-methylpyrrolidine,
    # is formed from the amine and cyclohexanecarbaldehyde. Vinylcyclohexane (an alkene) is incorrect.
    correct_reagent_A = "cyclohexanecarbaldehyde"

    # Principle 2: Enamine synthesis is an acid-catalyzed dehydration. A strong acid like TsOH
    # is a standard and highly effective catalyst for this type of reaction, especially with heat,
    # compared to a weak acid like acetic acid.
    correct_catalyst_B = "TsOH"

    # --- Step 2: Define the options as presented in the question ---
    options = {
        "A": {"A": "vinylcyclohexane", "B": "TsOH"},
        "B": {"A": "cyclohexanecarbaldehyde", "B": "TsOH"},
        "C": {"A": "vinylcyclohexane", "B": "Acetic acid"},
        "D": {"A": "cyclohexanecarbaldehyde", "B": "Acetic acid"}
    }

    # --- Step 3: Parse the LLM's final answer ---
    match = re.search(r'<<<([A-D])>>>', llm_answer)
    if not match:
        return "Failed to parse the answer. The answer should be in the format <<<X>>> where X is A, B, C, or D."

    selected_option_key = match.group(1)

    if selected_option_key not in options:
        return f"Invalid option '{selected_option_key}' selected. Valid options are A, B, C, D."

    selected_answer_components = options[selected_option_key]
    selected_reagent_A = selected_answer_components["A"]
    selected_catalyst_B = selected_answer_components["B"]

    # --- Step 4: Check if the selected answer satisfies the chemical constraints ---

    # Check Reagent A
    if selected_reagent_A != correct_reagent_A:
        return (f"Incorrect. The selected option '{selected_option_key}' is wrong because Reagent A is incorrect. "
                f"The reaction is an enamine synthesis, which requires a carbonyl compound like '{correct_reagent_A}', "
                f"not an alkene like '{selected_reagent_A}'.")

    # Check Catalyst B
    if selected_catalyst_B != correct_catalyst_B:
        return (f"Incorrect. The selected option '{selected_option_key}' is wrong because Catalyst B is incorrect. "
                f"While '{selected_catalyst_B}' is an acid, a strong acid like '{correct_catalyst_B}' is the standard and more "
                f"effective catalyst for this acid-catalyzed dehydration reaction, especially when heating is involved.")

    # --- Step 5: If all checks pass, the answer is correct ---
    return "Correct"

# Provided question and final answer from the prompt
question = """
A student performed an acid-catalyzed reaction of 3-methylpyrrolidine with compound A, the final product of the reaction is 1-(cyclohexylidenemethyl)-3-methylpyrrolidine. Select the suitable reagent (A) and catalyst (B) for the reaction.
3-methylpyrrolidine + A (B, Heat, solvent) ---> 1-(cyclohexylidenemethyl)-3-methylpyrrolidine

A) A = vinylcyclohexane, B = TsOH
B) A = cyclohexanecarbaldehyde, B = TsOH
C) A = vinylcyclohexane, B = Acetic acid
D) A = cyclohexanecarbaldehyde, B = Acetic acid
"""

final_answer = """
Here is a step-by-step analysis to determine the final answer based on the provided candidate responses.

**1. Analysis of the Chemical Reaction**

*   **Identify Reactants and Products:** The reaction starts with 3-methylpyrrolidine, which is a **secondary amine** (the nitrogen atom is bonded to two carbons and one hydrogen). The final product is 1-(cyclohexylidenemethyl)-3-methylpyrrolidine, which contains an `N-C=C` functional group, the defining feature of an **enamine**.
*   **Identify Reaction Type:** The transformation of a secondary amine into an enamine is a classic condensation reaction known as **enamine synthesis**. This reaction requires a carbonyl compound (an aldehyde or a ketone) and an acid catalyst.

**2. Determine the Identity of Reagent A**

*   Since the reaction is an enamine synthesis, Reagent A must be a carbonyl compound. The options provided are `vinylcyclohexane` (an alkene) and `cyclohexanecarbaldehyde` (an aldehyde).
*   An alkene is not the correct type of reactant. Therefore, Reagent A must be **cyclohexanecarbaldehyde**.
*   This is confirmed by a retrosynthetic analysis: The product `(pyrrolidine)-N-CH=C(cyclohexane)` is formed from the amine and a carbonyl compound. The carbon atom bonded to the nitrogen (`-CH=`) must have been the original carbonyl carbon (`C=O`). This points directly to cyclohexanecarbaldehyde.
*   **Conclusion for A:** Reagent A is cyclohexanecarbaldehyde. This eliminates options A and C.

**3. Determine the Identity of Catalyst B**

*   The reaction is an acid-catalyzed dehydration. The choice for the catalyst (B) is between `TsOH` (p-toluenesulfonic acid) and `Acetic acid`.
*   `TsOH` is a strong organic acid, while `Acetic acid` is a weak one.
*   The role of the acid is to facilitate the removal of a water molecule, which is often the slowest step. A strong acid like **TsOH** is a much more effective and standard catalyst for driving dehydration reactions to completion, especially when heat is applied (as indicated in the problem) to remove the water byproduct.
*   **Conclusion for B:** TsOH is the more suitable catalyst.

**4. Synthesize the Final Answer**

*   Based on the chemical principles, the correct combination is:
    *   Reagent A: **cyclohexanecarbaldehyde**
    *   Catalyst B: **TsOH**
*   Now, we map this correct combination to the options provided in the question:
    *   A) A = vinylcyclohexane, B = TsOH
    *   **B) A = cyclohexanecarbaldehyde, B = TsOH**
    *   C) A = vinylcyclohexane, B = Acetic acid
    *   D) A = cyclohexanecarbaldehyde, B = Acetic acid
*   The correct combination matches option B.

**5. Analysis of Candidate Answers**
A review of the 15 candidate answers reveals a crucial pattern:
*   **Reasoning Consensus:** Nearly every candidate correctly identifies that the reaction requires cyclohexanecarbaldehyde and that TsOH is the superior catalyst. The chemical logic is sound and consistent across the board.
*   **Conclusion Discrepancy:** The final letter choices (A, B, C, D) are scattered. This strongly suggests that the options may have been presented in a different order to the agents, causing them to map their correct reasoning to the wrong letter.

By trusting the overwhelming consensus on the chemical reasoning and mapping it to the original question's options, the correct answer is unequivocally determined.

<<<B>>>
"""

# Run the check
result = check_correctness(question, final_answer)
print(result)