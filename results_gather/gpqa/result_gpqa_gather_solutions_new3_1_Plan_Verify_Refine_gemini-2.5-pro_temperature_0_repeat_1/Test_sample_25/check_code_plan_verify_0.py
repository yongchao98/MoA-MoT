import re

def check_chemistry_answer(candidate_answer_text: str) -> str:
    """
    Checks the correctness of the answer for the given organic chemistry question.

    The function verifies two parts:
    1. The identity of reactant A.
    2. The reactivity order of the dienes (B).
    """

    # --- Part 1: Define the correct answer based on chemical principles ---

    # Part A: Identify Reactant A
    # The reaction is a [2+2] cycloaddition of cyclohexene with a ketene.
    # The product, 8,8-diiodobicyclo[4.2.0]octan-7-one, requires the ketene
    # to be diiodoketene (I2C=C=O).
    # The IUPAC name for I2C=C=O is 2,2-diiodoethen-1-one.
    correct_reactant_A = "2,2-diiodoethen-1-one"

    # Part B: Determine the reactivity order of dienes
    # Reactivity in Diels-Alder reactions depends on the ability to adopt the
    # s-cis conformation and electronic effects from substituents.
    # 3. cyclopenta-1,3-diene: Locked in s-cis conformation. Most reactive.
    # 1. 2,3-dimethylbuta-1,3-diene: Strong electron-donating groups (EDGs) on internal carbons, s-cis is accessible.
    # 2. (2E,4E)-hexa-2,4-diene: EDGs on terminal carbons, less activating than (1).
    # 4. (2Z,4Z)-hexa-2,4-diene: s-cis conformation is sterically impossible. Least reactive.
    # The correct order from most to least reactive is 3 > 1 > 2 > 4.
    correct_reactivity_B = [3, 1, 2, 4]

    # --- Part 2: Parse the candidate's answer ---

    # Find the final letter choice, e.g., <<<D>>>
    match = re.search(r'<<<([A-D])>>>', candidate_answer_text)
    if not match:
        return "Format Error: The answer does not contain a final choice in the format '<<<X>>>'."

    candidate_choice = match.group(1)

    # Define the options as provided in the question
    options = {
        'A': {"A": "4,4-diiodocyclobut-2-en-1-one", "B": [4, 2, 1, 3]},
        'B': {"A": "2,2-diiodoethen-1-one", "B": [4, 2, 1, 3]},
        'C': {"A": "4,4-diiodocyclobut-2-en-1-one", "B": [3, 1, 2, 4]},
        'D': {"A": "2,2-diiodoethen-1-one", "B": [3, 1, 2, 4]}
    }

    # Get the A and B values for the candidate's chosen option
    candidate_A = options[candidate_choice]["A"]
    candidate_B = options[candidate_choice]["B"]

    # --- Part 3: Check the correctness of the chosen option ---

    # Check Part A (Reactant)
    if candidate_A != correct_reactant_A:
        return (f"Incorrect. The identification of reactant A is wrong. "
                f"The answer chose '{candidate_A}', but the correct reactant is '{correct_reactant_A}'. "
                f"The reaction is a [2+2] cycloaddition between cyclohexene and diiodoketene.")

    # Check Part B (Reactivity Order)
    if candidate_B != correct_reactivity_B:
        return (f"Incorrect. The reactivity order of the dienes (B) is wrong. "
                f"The answer chose the order {candidate_B}, but the correct order is {correct_reactivity_B}. "
                f"The ranking (most to least reactive) is based on the ability to form the s-cis conformation and electronic effects: "
                f"cyclopenta-1,3-diene (3) > 2,3-dimethylbuta-1,3-diene (1) > (2E,4E)-hexa-2,4-diene (2) > (2Z,4Z)-hexa-2,4-diene (4).")

    # If both parts of the chosen option are correct, the answer is correct.
    return "Correct"

# The provided answer to be checked
llm_answer = """
An analysis of the provided information leads to the correct answer through a two-part logical deduction.

### **Step-by-step Analysis**

**Part 1: Identifying Reactant A**

1.  **Reaction Type:** The reaction starts with cyclohexene (a six-membered ring) and produces `8,8-diiodobicyclo[4.2.0]octan-7-one`. The product name `bicyclo[4.2.0]octan` signifies a structure where a six-membered ring is fused to a four-membered ring. This transformation from an alkene to a fused cyclobutane derivative is a classic example of a **[2+2] cycloaddition**.
2.  **Reactant Structure:** In this reaction, the two carbons of the cyclohexene double bond form two sides of the new four-membered ring. Reactant A must provide the other two carbons. The product's substituents, a ketone (`-7-one`) and two iodine atoms on the adjacent carbon (`8,8-diiodo`), must come from reactant A. For a [2+2] cycloaddition to occur, reactant A must contain a double bond. The only structure that fits these requirements is a ketene: `I₂C=C=O`.
3.  **Reactant Name:** The systematic IUPAC name for `I₂C=C=O` is **2,2-diiodoethen-1-one**.
4.  **Conclusion for Part A:** Based on this analysis, reactant A must be `2,2-diiodoethen-1-one`. This eliminates options A and C.

**Part 2: Determining the Reactivity Order of Dienes (B)**

The reactivity of dienes in cycloaddition reactions (like the Diels-Alder reaction) is governed by two primary factors:

*   **Conformation:** The diene must be in the planar *s-cis* conformation to react. Dienes that are locked in or can easily adopt this conformation are more reactive.
*   **Electronic Effects:** Electron-donating groups (EDGs) on the diene increase its reactivity.

Let's evaluate each diene based on these principles:

1.  **3. cyclopenta-1,3-diene:** This diene is cyclic, which rigidly locks its double bonds in the ideal *s-cis* conformation. This pre-organization for reaction makes it **the most reactive** diene in the list.
2.  **4. (2Z,4Z)-hexa-2,4-diene:** For this diene to adopt the required *s-cis* conformation, the two methyl groups would point inward, causing severe steric clash. This makes the reactive conformation energetically inaccessible. Therefore, this diene is essentially unreactive and is **the least reactive** in the list.
3.  **1. 2,3-dimethylbuta-1,3-diene:** This diene has two electron-donating methyl groups on its internal carbons (C2 and C3). These EDGs strongly activate the diene, and it can easily adopt the *s-cis* conformation. It is highly reactive.
4.  **2. (2E,4E)-hexa-2,4-diene:** This diene also has two activating methyl groups, but they are on the terminal carbons. While it can adopt the *s-cis* conformation, it is generally considered less reactive than diene (1) because internal substituents have a stronger activating effect.

**Reactivity Ranking (Most to Least Reactive):**
*   **Most Reactive:** 3 (cyclopenta-1,3-diene) - Locked *s-cis*.
*   **Second Most Reactive:** 1 (2,3-dimethylbuta-1,3-diene) - Strong electronic activation, easy *s-cis*.
*   **Third Most Reactive:** 2 ((2E,4E)-hexa-2,4-diene) - Less activation than (1).
*   **Least Reactive:** 4 ((2Z,4Z)-hexa-2,4-diene) - Sterically hindered from adopting *s-cis*.

The correct reactivity order is **3, 1, 2, 4**.

### **Final Conclusion**

*   Reactant **A** is **2,2-diiodoethen-1-one**.
*   The reactivity order **B** is **3, 1, 2, 4**.

Matching these two conclusions with the given options:
*   A) Incorrect A, Incorrect B
*   B) Correct A, Incorrect B
*   C) Incorrect A, Correct B
*   D) **Correct A, Correct B**

The correct option is D.

<<<D>>>
"""

# print(check_chemistry_answer(llm_answer))