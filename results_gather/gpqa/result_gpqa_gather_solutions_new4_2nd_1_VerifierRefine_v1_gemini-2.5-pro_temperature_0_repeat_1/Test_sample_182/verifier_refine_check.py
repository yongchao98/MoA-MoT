import re

def check_chemistry_answer(llm_answer_text: str) -> str:
    """
    Checks the correctness of the LLM's answer for the given organic chemistry question.

    The question asks for the Index of Hydrogen Deficiency (IHD) of the product obtained
    when 2-formyl-5-vinylcyclohex-3-enecarboxylic acid is reacted with red phosphorus and excess HI.
    """
    try:
        # --- Step 1: Define the chemical principles and problem parameters ---

        # The starting material is 2-formyl-5-vinylcyclohex-3-enecarboxylic acid.
        # Let's determine the sources of unsaturation (IHD = rings + pi bonds).
        # - `cyclohex-`: 1 ring
        # - `-3-ene`: 1 pi bond (C=C)
        # - `vinyl`: 1 pi bond (C=C)
        # - `formyl`: 1 pi bond (C=O)
        # - `carboxylic acid`: 1 pi bond (C=O)
        
        # The reaction is with red phosphorus and excess HI. This is a powerful reduction.
        # It reduces all pi bonds (C=C and C=O) to single bonds.
        # It does not break the stable ring structure.

        # --- Step 2: Determine the structure and IHD of the product ---

        # The number of rings in the product remains the same as in the starting material.
        num_rings_product = 1

        # All pi bonds are reduced by the reaction.
        num_pi_bonds_product = 0

        # Calculate the correct IHD for the product.
        correct_ihd = num_rings_product + num_pi_bonds_product

        # --- Step 3: Parse the LLM's final answer from the text ---
        match = re.search(r'<<<([A-D])>>>', llm_answer_text)
        if not match:
            return "Failure: Could not find the final answer in the required format '<<<X>>>'."

        llm_choice = match.group(1)

        # --- Step 4: Map the LLM's choice to its numerical value based on the question's options ---
        options = {'A': 1, 'B': 0, 'C': 3, 'D': 5}
        if llm_choice not in options:
            return f"Failure: The selected option '{llm_choice}' is not a valid choice (A, B, C, or D)."

        llm_ihd_value = options[llm_choice]

        # --- Step 5: Compare the LLM's answer with the calculated correct answer ---
        if llm_ihd_value == correct_ihd:
            return "Correct"
        else:
            reason = (
                f"The answer is incorrect.\n"
                f"Reasoning:\n"
                f"1. The starting material is 2-formyl-5-vinylcyclohex-3-enecarboxylic acid.\n"
                f"2. The reaction with red phosphorus and excess HI is a complete reduction. It reduces all pi bonds (C=C and C=O) but does not break the ring.\n"
                f"3. The product is a saturated substituted cyclohexane. It has 1 ring and 0 pi bonds.\n"
                f"4. The Index of Hydrogen Deficiency (IHD) of the product is the sum of rings and pi bonds, which is 1 + 0 = 1.\n"
                f"5. The provided answer '{llm_choice}' corresponds to an IHD of {llm_ihd_value}, but the correct IHD is {correct_ihd}."
            )
            return reason

    except Exception as e:
        return f"An error occurred during the check: {e}"

# The user-provided text containing the LLM's answer.
llm_response = """
Here is a step-by-step analysis of the question to determine the final answer.

**1. Deconstruct the Starting Material**
The starting material is **2-formyl-5-vinylcyclohex-3-enecarboxylic acid**. To understand the reaction, we must identify all its structural features, particularly those contributing to the Index of Hydrogen Deficiency (IHD). The IHD is the sum of the number of rings and pi (π) bonds in a molecule.

*   **Ring Structure:** The molecule has a six-membered ring (`cyclohex-`). This accounts for **1** degree of unsaturation.
*   **Pi Bonds (π bonds):**
    *   One carbon-carbon double bond (C=C) is present within the ring (`-3-ene`). This adds **1** to the IHD.
    *   One carbon-carbon double bond (C=C) is in the vinyl group (`-CH=CH₂`). This adds **1** to the IHD.
    *   One carbon-oxygen double bond (C=O) is in the formyl group (`-CHO`). This adds **1** to the IHD.
    *   One carbon-oxygen double bond (C=O) is in the carboxylic acid group (`-COOH`). This adds **1** to the IHD.

The total IHD of the starting material is 1 (ring) + 4 (π bonds) = 5. This is a potential distractor, as the question asks for the IHD of the *product*.

**2. Analyze the Reaction and Reagents**
The reaction is carried out with **red phosphorus and an excess of HI (hydroiodic acid)**. This combination is one of the most powerful reducing agents in organic chemistry. Its primary function is to perform a complete reduction of various functional groups to their corresponding saturated hydrocarbon (alkane) forms.

*   **Effect on Functional Groups:**
    *   It reduces aldehydes (-CHO) and carboxylic acids (-COOH) completely to methyl groups (-CH₃).
    *   It reduces carbon-carbon double bonds (C=C) to single bonds (C-C).
*   **Effect on the Carbon Skeleton:** The reaction is not strong enough to break the stable cycloalkane ring. The ring structure remains intact.

**3. Determine the Structure of the Final Product**
The reaction will reduce all the pi bonds identified in Step 1 but will leave the ring structure itself untouched.

*   The C=C bond in the ring is reduced to a C-C single bond.
*   The C=C bond in the vinyl group is reduced, converting it to an ethyl group (-CH₂CH₃).
*   The C=O bond in the formyl group is eliminated, reducing the group to a methyl group (-CH₃).
*   The C=O bond in the carboxylic acid group is eliminated, reducing the group to a methyl group (-CH₃).

The final product is a **fully saturated, substituted cyclohexane**. It contains only single bonds and one ring.

**4. Calculate the Index of Hydrogen Deficiency (IHD) of the Product**
The IHD of the product is the sum of its rings and pi bonds.

*   **Number of Rings:** The cyclohexane ring is preserved. The product has **1** ring.
*   **Number of Pi Bonds:** All double bonds (both C=C and C=O) were reduced to single bonds. The product has **0** pi bonds.

Therefore, the IHD of the final product is:
IHD = (Number of Rings) + (Number of Pi Bonds) = 1 + 0 = **1**.

This single degree of unsaturation is due to the presence of the ring.

**5. Conclusion**
The analysis consistently shows that the reaction reduces all pi bonds but preserves the single ring structure. The Index of Hydrogen Deficiency of the final product is 1. This corresponds to option A.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(llm_response)
print(result)