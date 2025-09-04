import re

def check_chemistry_answer(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given chemistry question.
    The function verifies both the final selected option and the logical reasoning provided.
    """

    # --- Define Ground Truth based on chemical principles ---
    # Reaction: Stork enamine synthesis from a secondary amine and an aldehyde.
    # Reagent A (carbonyl): Must be cyclohexanecarbaldehyde to form the cyclohexylidenemethyl group.
    # Reagent B (catalyst): Must be an acid. Strong acids like TsOH are preferred for efficient dehydration.
    correct_option = 'C'
    correct_reagent_A = "cyclohexanecarbaldehyde"
    incorrect_reagent_A = "vinylcyclohexane"
    correct_catalyst_B = "TsOH"
    incorrect_catalyst_B = "Acetic acid"

    # --- Check 1: Final Answer Extraction and Comparison ---
    # The final answer should be enclosed in <<< >>>
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The answer is not in the required format '<<<...>>>'."
    
    llm_choice = match.group(1)
    if llm_choice != correct_option:
        return f"Incorrect final answer. The selected option is {llm_choice}, but the correct option is {correct_option} based on chemical principles."

    # --- Check 2: Verification of Reasoning ---
    # The reasoning should correctly identify the key components and principles.
    
    # 2a. Check for correct identification of the reaction type (enamine synthesis).
    if "enamine" not in llm_answer_text.lower():
        return "Incorrect reasoning: The explanation fails to identify the product as an enamine or the reaction as an enamine synthesis."

    # 2b. Check for correct identification of Reagent A.
    if correct_reagent_A not in llm_answer_text:
        return f"Incorrect reasoning: The explanation fails to identify the correct Reagent A, which is {correct_reagent_A}."
    
    # 2c. Check for correct elimination of the incorrect Reagent A.
    if incorrect_reagent_A in llm_answer_text and "eliminate" not in llm_answer_text.lower():
        return f"Incorrect reasoning: The explanation mentions {incorrect_reagent_A} but does not correctly state that it should be eliminated."

    # 2d. Check for correct identification of Catalyst B and the justification.
    if correct_catalyst_B not in llm_answer_text:
        return f"Incorrect reasoning: The explanation fails to identify the correct Catalyst B, which is {correct_catalyst_B}."
    
    # 2e. Check for the justification of catalyst choice (strong vs. weak acid).
    if not ("strong acid" in llm_answer_text.lower() and "tsoh" in llm_answer_text.lower()):
        return f"Incorrect reasoning: The explanation fails to justify the choice of {correct_catalyst_B} by identifying it as a strong acid suitable for the reaction."
    
    if not ("weak acid" in llm_answer_text.lower() and "acetic acid" in llm_answer_text.lower()):
        return f"Incorrect reasoning: The explanation fails to provide a comparison by identifying {incorrect_catalyst_B} as a weaker, less suitable catalyst."

    # If all checks pass, the answer and its reasoning are correct.
    return "Correct"

# The user-provided LLM answer
llm_answer = """
Here is a step-by-step analysis of the problem:

**Step 1: Analyze the Product Structure**

The product is named 1-(cyclohexylidenemethyl)-3-methylpyrrolidine. Let's break down this name to understand its structure.
*   **3-methylpyrrolidine:** This is the parent amine structure, a five-membered ring containing nitrogen, with a methyl group at the 3-position. Since it's a secondary amine, the nitrogen atom has one hydrogen attached.
*   **1-(...):** This indicates that a substituent is attached to the nitrogen atom (position 1) of the pyrrolidine ring.
*   **cyclohexylidenemethyl:** This is the substituent. It consists of a "cyclohexylidene" group and a "methyl" group.
    *   A "cyclohexylidene" group is a cyclohexane ring attached to another atom via a double bond (`=C6H10`).
    *   A "methyl" group in this context refers to a `-CH=` fragment.
    *   Combining them, the `cyclohexylidenemethyl` group is `-CH=C6H10`.
*   **Putting it all together:** The product has the 3-methylpyrrolidine group attached via its nitrogen atom to the `-CH=` carbon of the `cyclohexylidenemethyl` group. This forms a `N-C=C` linkage, which is the defining feature of an **enamine**.

**Step 2: Identify the Reaction Type**

The reaction forms an enamine from a secondary amine (3-methylpyrrolidine). This is a classic organic reaction known as the **Stork enamine synthesis**. It is a condensation reaction between a secondary amine and a carbonyl compound (an aldehyde or a ketone). The reaction is acid-catalyzed and involves the removal of a water molecule.

**Step 3: Determine the Starting Carbonyl Compound (Reagent A)**

To find the carbonyl compound (A) that reacted with 3-methylpyrrolidine, we can perform a "retrosynthesis" by conceptually hydrolyzing the enamine product.
*   The enamine structure is: `(3-methylpyrrolidinyl)-N-CH=C(cyclohexane ring)`
*   The reaction forms the C=C double bond by removing a proton from the carbon alpha (adjacent) to the carbonyl group of the starting material. The C=O group itself becomes the C-N bond.
*   Let's reverse the process:
    1.  Cleave the N-C bond of the enamine.
    2.  Add a hydrogen atom back to the nitrogen to reform the secondary amine (3-methylpyrrolidine).
    3.  Place a carbonyl group (C=O) on the carbon that was attached to the nitrogen.
    4.  Convert the C=C double bond back to a C-C single bond by adding a hydrogen to the other carbon of the original double bond.
*   Applying this to our product:
    *   The carbon attached to the nitrogen is the `CH` group. This becomes the carbonyl carbon, `C=O`.
    *   The other carbon of the double bond is part of the cyclohexane ring. It gets a hydrogen atom back.
    *   So, the starting carbonyl compound is `O=CH-(cyclohexane ring)`.
*   This structure is **cyclohexanecarbaldehyde**.
*   Looking at the options for A, we can eliminate A) and B) which list vinylcyclohexane. The correct reagent A is cyclohexanecarbaldehyde. This leaves us with options C and D.

**Step 4: Select the Suitable Catalyst (Reagent B)**

The reaction is an acid-catalyzed dehydration. The role of the acid is to protonate the hydroxyl group of the carbinolamine intermediate, turning it into a good leaving group (water).
*   The choices for the catalyst (B) are **TsOH** (p-toluenesulfonic acid) and **Acetic acid**.
*   **TsOH** is a strong organic acid (pKa ≈ -2.8).
*   **Acetic acid** is a weak organic acid (pKa ≈ 4.76).
*   For enamine formation, a catalytic amount of a strong acid like TsOH is generally preferred. It is effective at promoting the dehydration step without being so concentrated that it fully protonates the starting amine (which would make it non-nucleophilic) or the enamine product (which would reverse the reaction). TsOH is the standard catalyst for this type of reaction.
*   While acetic acid can catalyze the reaction, it is much weaker and therefore less efficient. The dehydration step would be significantly slower.
*   Therefore, **TsOH** is the more suitable and commonly used catalyst for this transformation.

**Step 5: Final Conclusion**

*   Reagent A must be **cyclohexanecarbaldehyde**.
*   Catalyst B should be a strong acid catalyst, making **TsOH** the better choice.
*   Combining these findings, the correct option is C.

<<<C>>>
"""

# Run the check
result = check_chemistry_answer(llm_answer)
print(result)