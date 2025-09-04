import re

def check_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for a multi-step organic chemistry problem.

    The function verifies two main aspects:
    1. The chemical reasoning to determine the final product and its properties.
    2. The final answer choice based on the derived properties.

    Args:
        llm_answer_text: A string containing the full text of the LLM's answer,
                         which should include the final choice in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or an explanation of why the answer is incorrect.
    """

    # Step 1: Determine the correct answer based on chemical principles.
    # The reaction sequence is:
    # 1. Cyclohexanone + Br2 -> 2-bromocyclohexanone (alpha-bromination)
    # 2. 2-bromocyclohexanone + NaOH -> Cyclopentanecarboxylic acid (Favorskii rearrangement)
    # 3. Cyclopentanecarboxylic acid + SOCl2 -> Cyclopentanecarbonyl chloride (acid chloride formation)
    # 4. Cyclopentanecarbonyl chloride + LiAlH(O-t-Bu)3 -> Cyclopentanecarbaldehyde (selective reduction)
    # The final product is Cyclopentanecarbaldehyde.

    # Counting chemically distinct hydrogens in Cyclopentanecarbaldehyde:
    # The molecule has a plane of symmetry passing through the C1-CHO bond.
    # - Aldehyde proton: 1 type
    # - C1 proton (on the ring): 1 type
    # - C2/C5 protons: These 4 protons are diastereotopic, forming 2 distinct types.
    # - C3/C4 protons: These 4 protons are also diastereotopic, forming another 2 distinct types.
    # Total = 1 + 1 + 2 + 2 = 6
    correct_numerical_answer = 6

    # Step 2: Define the mapping from options to numbers as specified in the prompt's analysis.
    # Question options: A) 10, B) 7, C) 6, D) 8
    options_map = {
        'A': 10,
        'B': 7,
        'C': 6,
        'D': 8
    }

    # Step 3: Extract the chosen option from the LLM's answer text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect: The answer is not in the required format '<<<X>>>' where X is one of A, B, C, or D."

    chosen_option = match.group(1)
    
    # Step 4: Get the numerical value corresponding to the chosen option.
    llm_numerical_answer = options_map.get(chosen_option)
    if llm_numerical_answer is None:
        return f"Incorrect: The chosen option '{chosen_option}' is invalid."

    # Step 5: Compare the LLM's answer with the correct answer and provide a reason if incorrect.
    if llm_numerical_answer == correct_numerical_answer:
        # The reasoning in the text also correctly identifies the product and counts 6 hydrogens.
        # The final choice <<<C>>> correctly maps to 6.
        return "Correct"
    else:
        reason = (
            f"Incorrect: The final answer is {chosen_option}, which corresponds to {llm_numerical_answer} distinct hydrogens according to the provided option mapping.\n"
            f"The correct number of chemically distinct hydrogen atoms in the final product (cyclopentanecarbaldehyde) is {correct_numerical_answer}.\n"
            "The correct analysis is as follows:\n"
            "- 1 type for the aldehyde proton.\n"
            "- 1 type for the proton on the ring carbon attached to the aldehyde group (C1).\n"
            "- 2 types for the four diastereotopic protons on the adjacent carbons (C2 and C5).\n"
            "- 2 types for the four diastereotopic protons on the next carbons (C3 and C4).\n"
            f"The sum is 1 + 1 + 2 + 2 = {correct_numerical_answer}, which corresponds to option C."
        )
        return reason

# The user-provided answer to be checked
llm_answer = """
Based on a step-by-step analysis of the chemical problem and a review of the candidate answers, a definitive conclusion can be reached.

### Step-by-Step Analysis

**1. Determine the Final Product (Product 4)**

*   **Step 1: Cyclohexanone + Bromine → Product 1**
    This is a standard α-bromination of a ketone, yielding **2-bromocyclohexanone**.
*   **Step 2: Product 1 + NaOH, heat → Product 2**
    This is the reaction of an α-halo ketone with a strong base. The conditions are classic for a **Favorskii rearrangement**, which involves a ring contraction. The six-membered ring rearranges to form a five-membered ring carboxylic acid. The product is **cyclopentanecarboxylic acid**.
*   **Step 3: Product 2 + Thionyl chloride (SOCl₂) → Product 3**
    This is a standard conversion of a carboxylic acid to a more reactive acyl chloride. The product is **cyclopentanecarbonyl chloride**. This step confirms that the Favorskii rearrangement was the correct pathway in Step 2.
*   **Step 4: Product 3 + Lithium tri-tert-butoxyaluminum hydride → Product 4**
    This reagent is a mild, sterically hindered reducing agent specifically used to reduce acyl chlorides to aldehydes without over-reduction to an alcohol. The final product is **cyclopentanecarbaldehyde**.

There is a unanimous consensus among all candidate answers that the final product is cyclopentanecarbaldehyde.

**2. Count the Chemically Distinct Hydrogen Atoms in Product 4**

The task is to analyze the structure of **cyclopentanecarbaldehyde** and count its unique hydrogen environments, which corresponds to the number of signals in a ¹H NMR spectrum.

*   **Structure and Symmetry:** The molecule has a plane of symmetry that passes through the aldehyde group (-CHO), the ring carbon it is attached to (C1), and the hydrogen on C1. This plane bisects the C3-C4 bond.
*   **Counting the Unique Hydrogen Environments:**
    1.  **The Aldehyde Proton:** The single hydrogen of the -CHO group is in a unique environment. (1 type)
    2.  **The C1 Proton:** The single hydrogen on the ring carbon attached to the aldehyde group is also unique. (1 type)
    3.  **The C2 and C5 Protons:** The plane of symmetry makes the CH₂ groups at C2 and C5 equivalent. However, the two protons on each of these carbons are **diastereotopic** because of their different spatial relationships to the aldehyde group on the adjacent prochiral carbon (C1). This results in two distinct types of hydrogens from these four protons. (2 types)
    4.  **The C3 and C4 Protons:** Similarly, the CH₂ groups at C3 and C4 are equivalent. The two protons on each of these carbons are also diastereotopic, resulting in another two distinct types of hydrogens. (2 types)

*   **Total Count:** Summing the distinct types gives: 1 + 1 + 2 + 2 = **6**.

### Careful Points and Synthesis

1.  **Key Reaction:** The correct identification of the **Favorskii rearrangement** in Step 2 is the most critical part of determining the reaction pathway. An alternative elimination reaction is inconsistent with the subsequent steps.
2.  **Reagent Specificity:** Recognizing that Lithium tri-tert-butoxyaluminum hydride is a selective reducing agent that stops at the aldehyde is crucial. A stronger reagent like LiAlH₄ would produce an alcohol, leading to a different answer.
3.  **Symmetry and Diastereotopicity:** The final count of 6 relies on correctly identifying the molecule's plane of symmetry and understanding the concept of diastereotopic protons. Even with a plane of symmetry, protons on a methylene group (CH₂) can be non-equivalent if adjacent to a prochiral center.
4.  **Mapping to Options:** The calculated number of distinct hydrogens is 6. According to the question's options (A=10, B=7, C=6, D=8), this corresponds to option C. Many of the candidate answers made clerical errors in their final letter selection despite performing the chemical analysis correctly.

The consistent and correct chemical analysis across all detailed answers confirms that there are 6 chemically distinct hydrogen atoms in the final product.

<<<C>>>
"""

# print(check_chemistry_answer(llm_answer))