import re

def check_correctness(question, llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The function verifies the choice of reagent (A) and catalyst (B) based on
    the principles of enamine synthesis.

    Args:
        question (str): The question text.
        llm_answer_text (str): The text of the LLM's answer, including the final choice.

    Returns:
        str: "Correct" if the answer is correct, otherwise a reason for the error.
    """

    # --- Define the correct chemical principles for this reaction ---
    # 1. Reaction type: Enamine synthesis from a secondary amine.
    # 2. Reagent A must be a carbonyl compound that yields the product structure.
    #    Retrosynthesis of 1-(cyclohexylidenemethyl)-3-methylpyrrolidine points to
    #    cyclohexanecarbaldehyde. Vinylcyclohexane is an alkene and is incorrect.
    correct_reagent_A = "cyclohexanecarbaldehyde"

    # 3. Catalyst B must be a suitable acid for dehydration.
    #    TsOH is a strong acid and a standard, highly effective catalyst for
    #    dehydration reactions, especially when heat is applied. Acetic acid is
    #    a weak acid and is far less effective, making it a less suitable choice.
    correct_catalyst_B = "TsOH"

    # --- Define the options from the question ---
    options = {
        "A": {"reagent_A": "vinylcyclohexane", "catalyst_B": "Acetic acid"},
        "B": {"reagent_A": "vinylcyclohexane", "catalyst_B": "TsOH"},
        "C": {"reagent_A": "cyclohexanecarbaldehyde", "catalyst_B": "TsOH"},
        "D": {"reagent_A": "cyclohexanecarbaldehyde", "catalyst_B": "Acetic acid"}
    }

    # --- Extract the final answer from the LLM's response ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find a final answer in the format <<<A>>>, <<<B>>>, etc."

    final_answer_key = match.group(1)
    
    if final_answer_key not in options:
        return f"Error: The extracted answer '{final_answer_key}' is not a valid option."

    selected_option = options[final_answer_key]

    # --- Check the selected option against the correct principles ---
    # Check Reagent A
    if selected_option["reagent_A"] != correct_reagent_A:
        return (f"Incorrect. The selected answer is {final_answer_key}, which proposes "
                f"Reagent A = {selected_option['reagent_A']}. "
                f"The reaction is an enamine synthesis, which requires a carbonyl compound. "
                f"The correct Reagent A is '{correct_reagent_A}'.")

    # Check Catalyst B
    if selected_option["catalyst_B"] != correct_catalyst_B:
        return (f"Incorrect. The selected answer is {final_answer_key}, which proposes "
                f"Catalyst B = {selected_option['catalyst_B']}. "
                f"For an efficient acid-catalyzed dehydration reaction, a strong acid like '{correct_catalyst_B}' "
                f"is the standard and most suitable catalyst, especially when heat is used. "
                f"'{selected_option['catalyst_B']}' is a less effective weak acid.")

    # If all checks pass, the answer is correct
    return "Correct"

# Provided question and LLM answer
question = """
A student performed an acid-catalyzed reaction of 3-methylpyrrolidine with compound A, the final product of the reaction is 1-(cyclohexylidenemethyl)-3-methylpyrrolidine. Select the suitable reagent (A) and catalyst (B) for the reaction.
3-methylpyrrolidine + A (B, Heat, solvent) ---> 1-(cyclohexylidenemethyl)-3-methylpyrrolidine

A) A = vinylcyclohexane, B = Acetic acid
B) A = vinylcyclohexane, B = TsOH
C) A = cyclohexanecarbaldehyde, B = TsOH
D) A = cyclohexanecarbaldehyde, B = Acetic acid
"""

llm_answer = """
Here is a thinking process to arrive at the final answer.

*   **Plan:**
    1.  First, I will analyze the reasoning within all the provided candidate answers to establish a consensus on the identity of Reagent A.
    2.  Next, I will analyze the reasoning for Catalyst B to determine the most suitable choice based on the chemical principles described.
    3.  Finally, I will combine these two findings to select the correct option and provide the final answer.

I will now execute the first step of the plan.

*   **Execution:**
    I have analyzed the reasoning for Reagent A across all 17 candidate answers.
    *   **Reaction Type:** All answers correctly identify the reaction as an enamine synthesis, which involves the condensation of a secondary amine with a carbonyl compound (an aldehyde or a ketone).
    *   **Reagent A Identity:** All answers correctly conclude that to form the specified enamine product, the secondary amine must react with a carbonyl compound. They correctly identify **cyclohexanecarbaldehyde** as the necessary aldehyde. They all correctly rule out vinylcyclohexane, as it is an alkene and would not participate in this reaction.
    *   **Conclusion:** There is a unanimous consensus in the *reasoning* of all 17 answers that Reagent A must be **cyclohexanecarbaldehyde**. This eliminates options A and D.

I will now refine the plan and proceed to the next step.

*   **Plan:**
    1.  Reagent A has been identified as cyclohexanecarbaldehyde.
    2.  Now, I will analyze the reasoning for Catalyst B, choosing between TsOH and Acetic acid (options C and D).
    3.  Finally, I will combine the findings to select the correct option.

I will now execute the second step of the plan.

*   **Execution:**
    I have analyzed the reasoning for Catalyst B across all 17 candidate answers.
    *   **Catalyst Role:** All answers correctly state that the acid catalyst's role is to facilitate the dehydration (elimination of water) step, which is often rate-limiting.
    *   **Catalyst Comparison:** The choice is between TsOH (a strong acid) and Acetic acid (a weak acid). The reasoning in every single candidate answer concludes that **TsOH** is the more suitable, effective, and standard catalyst for this type of dehydration reaction. It is described as a powerful catalyst that efficiently drives the reaction to completion, which is consistent with the use of heat. Acetic acid is consistently described as a weaker and less effective alternative.
    *   **Conclusion:** There is a unanimous consensus in the *reasoning* of all 17 answers that Catalyst B should be **TsOH**.

I will now refine the plan and proceed to the final step.

*   **Plan:**
    1.  Reagent A is cyclohexanecarbaldehyde.
    2.  Catalyst B is TsOH.
    3.  I will now combine these conclusions to select the final answer.

I am now ready to provide the final answer.

*   **Execution:**
    Based on the step-by-step analysis of the reasoning provided in all candidate answers:
    *   The correct Reagent A is **cyclohexanecarbaldehyde**.
    *   The most suitable Catalyst B is **TsOH**.
    *   This combination corresponds to option C.

<<<C>>>
"""

# Run the check
result = check_correctness(question, llm_answer)
print(result)