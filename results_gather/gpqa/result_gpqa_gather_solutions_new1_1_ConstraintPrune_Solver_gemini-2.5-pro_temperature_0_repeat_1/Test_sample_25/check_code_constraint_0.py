import re

def check_correctness_of_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The function verifies two parts of the answer:
    1. The identity of reactant A.
    2. The reactivity order of the dienes (B).

    Args:
        llm_answer_text: The full text of the LLM's response, including the final answer
                         in the format <<<X>>>.

    Returns:
        A string indicating "Correct" or explaining the reason for the incorrectness.
    """

    # --- Define the ground truth based on established chemical principles ---

    # Part A: Identity of Reactant A
    # The reaction is a [2+2] cycloaddition between cyclohexene and a ketene.
    # The product, 8,8-diiodobicyclo[4.2.0]octan-7-one, requires the ketene to be
    # diiodoketene (I2C=C=O). Its systematic IUPAC name is 2,2-diiodoethen-1-one.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # Part B: Reactivity Order of Dienes
    # The reactivity in Diels-Alder reactions is governed by the ability to form an s-cis
    # conformation and by electronic effects.
    # 1. cyclopenta-1,3-diene (3): Locked in the reactive s-cis conformation. Most reactive.
    # 2. 2,3-dimethylbuta-1,3-diene (1): Has activating electron-donating groups (EDGs) on
    #    internal carbons and can easily adopt the s-cis form. Very reactive.
    # 3. (2E,4E)-hexa-2,4-diene (2): Has terminal EDGs and can adopt the s-cis form, but is
    #    less reactive than (1).
    # 4. (2Z,4Z)-hexa-2,4-diene (4): Severe steric hindrance between the methyl groups
    #    prevents the formation of the s-cis conformation. Least reactive.
    correct_reactivity_order_B = [3, 1, 2, 4]

    # --- Define the options provided in the question ---
    options = {
        'A': {
            "reactant_A": "2,2-diiodoethen-1-one",
            "reactivity_B": [4, 2, 1, 3]
        },
        'B': {
            "reactant_A": "4,4-diiodocyclobut-2-en-1-one",
            "reactivity_B": [3, 1, 2, 4]
        },
        'C': {
            "reactant_A": "2,2-diiodoethen-1-one",
            "reactivity_B": [3, 1, 2, 4]
        },
        'D': {
            "reactant_A": "4,4-diiodocyclobut-2-en-1-one",
            "reactivity_B": [4, 2, 1, 3]
        }
    }

    # --- Extract the proposed answer from the LLM's text ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Error: Could not find the final answer in the required format '<<<X>>>' in the provided text."

    proposed_answer_key = match.group(1)
    if proposed_answer_key not in options:
        return f"Error: The extracted answer key '{proposed_answer_key}' is not a valid option (A, B, C, or D)."

    proposed_answer = options[proposed_answer_key]

    # --- Compare the proposed answer with the ground truth ---
    is_A_correct = (proposed_answer["reactant_A"] == correct_reactant_A)
    is_B_correct = (proposed_answer["reactivity_B"] == correct_reactivity_order_B)

    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = []
        if not is_A_correct:
            error_messages.append(
                f"Constraint for reactant A is not satisfied. "
                f"The proposed reactant A is '{proposed_answer['reactant_A']}', but the correct reactant is '{correct_reactant_A}'. "
                "The reaction is a [2+2] cycloaddition between cyclohexene and diiodoketene (I2C=C=O) to form the bicyclo[4.2.0]octanone product."
            )
        if not is_B_correct:
            error_messages.append(
                f"Constraint for diene reactivity B is not satisfied. "
                f"The proposed reactivity order is {proposed_answer['reactivity_B']}, but the correct order is {correct_reactivity_order_B}. "
                "The correct order (most to least reactive) is: cyclopenta-1,3-diene (3) > 2,3-dimethylbuta-1,3-diene (1) > (2E,4E)-hexa-2,4-diene (2) > (2Z,4Z)-hexa-2,4-diene (4)."
            )
        return "\n".join(error_messages)

# The user-provided text containing the final answer to be checked.
llm_response = """
Based on the provided candidate answers, here is a step-by-step analysis to determine the final correct answer.

The problem consists of two independent parts:
1.  **Part A:** Identifying the reactant 'A' in a specific chemical reaction.
2.  **Part B:** Determining the correct reactivity order for a list of four dienes.

The final answer must be correct for both parts.

### **Analysis of Part A: Identifying Reactant A**

1.  **Reaction Analysis:** The given reaction is `Cyclohexene + A ---> 8,8-diiodobicyclo[4.2.0]octan-7-one`.
2.  **Product Structure:** The product name `bicyclo[4.2.0]octan` signifies a fused ring system composed of a six-membered ring and a four-membered ring. Since the starting material is cyclohexene (a six-membered ring), the reaction must involve the formation of a new four-membered ring onto the cyclohexene backbone.
3.  **Reaction Type:** The formation of a four-membered ring from two double-bonded species is a **[2+2] cycloaddition**. The double bond of cyclohexene reacts with a double bond in reactant A.
4.  **Deducing Reactant A:** The product has a ketone (`-7-one`) and two iodine atoms (`8,8-diiodo`) on the newly formed four-membered ring. This means reactant A must provide a two-carbon unit containing these functional groups. For a [2+2] cycloaddition, this unit must contain a C=C double bond. The class of compounds that fits this description and is well-known to undergo [2+2] cycloadditions with alkenes is **ketenes** (R₂C=C=O). The specific ketene required here is diiodoketene, with the structure `I₂C=C=O`.
5.  **Nomenclature:** The systematic IUPAC name for `I₂C=C=O` is **2,2-diiodoethen-1-one**.
6.  **Conclusion for Part A:** Reactant A is `2,2-diiodoethen-1-one`. This eliminates options B and D, leaving A and C as possibilities.

### **Analysis of Part B: Diene Reactivity Order**

1.  **Governing Principles:** The reactivity of dienes in cycloaddition reactions (like the Diels-Alder reaction) is primarily determined by two factors:
    *   **Conformation:** The diene must be in the planar *s-cis* conformation to react. Dienes that are locked in this conformation are highly reactive, while those that are sterically hindered from achieving it are unreactive.
    *   **Electronic Effects:** Electron-donating groups (EDGs) on the diene increase its reactivity.
2.  **Evaluating the Dienes:**
    *   **3. cyclopenta-1,3-diene:** This diene is cyclic and is permanently locked in the required *s-cis* conformation. This makes it exceptionally reactive. It is the **most reactive**.
    *   **4. (2Z,4Z)-hexa-2,4-diene:** For this diene to adopt the *s-cis* conformation, the two terminal methyl groups would point inward and clash severely. This extreme steric hindrance makes the reactive conformation energetically inaccessible. It is the **least reactive**.
    *   **1. 2,3-dimethylbuta-1,3-diene:** This diene has two activating EDGs (methyl groups) on its internal carbons. It can easily adopt the *s-cis* conformation. It is highly reactive.
    *   **2. (2E,4E)-hexa-2,4-diene:** This diene has activating EDGs on its terminal carbons. It can adopt the *s-cis* conformation, but it is less favored than for diene (1) due to minor steric interactions and less effective electronic activation compared to internal substitution.
3.  **Reactivity Order:** Based on the analysis, the order from most reactive to least reactive is: `cyclopenta-1,3-diene` (3) > `2,3-dimethylbuta-1,3-diene` (1) > `(2E,4E)-hexa-2,4-diene` (2) > `(2Z,4Z)-hexa-2,4-diene` (4).
4.  **Conclusion for Part B:** The correct reactivity sequence is **3, 1, 2, 4**.

### **Final Conclusion**

*   From Part A, reactant A is **2,2-diiodoethen-1-one**.
*   From Part B, the reactivity order is **3, 1, 2, 4**.

We must find the option that matches both conclusions.
*   A) A = 2,2-diiodoethen-1-one, B = 4, 2, 1, 3 (Incorrect B)
*   B) A = 4,4-diiodocyclobut-2-en-1-one, B = 3, 1, 2, 4 (Incorrect A)
*   **C) A = 2,2-diiodoethen-1-one, B = 3, 1, 2, 4 (Correct A and Correct B)**
*   D) A = 4,4-diiodocyclobut-2-en-1-one, B = 4, 2, 1, 3 (Incorrect A and B)

The only option that satisfies both constraints is C. A careful review of the candidate answers shows that while most correctly deduced the chemical principles, many made a final error in matching their correct reasoning to the letter of the multiple-choice option.

<<<C>>>
"""

# Run the check
result = check_correctness_of_chemistry_answer(llm_response)
print(result)