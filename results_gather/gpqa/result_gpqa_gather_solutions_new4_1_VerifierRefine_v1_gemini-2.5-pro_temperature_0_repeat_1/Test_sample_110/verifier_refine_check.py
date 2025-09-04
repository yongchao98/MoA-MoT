import re

def check_answer_correctness(llm_answer_text: str) -> str:
    """
    Checks the correctness of the selected answer for the given chemistry question.

    Args:
        llm_answer_text: The full text of the LLM's response, including the final answer in <<<X>>> format.

    Returns:
        "Correct" if the answer is correct, otherwise a string explaining the reason for the error.
    """
    # Step 1: Define the correct products based on chemical principles.
    # Reaction A: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK)
    # This is a Michael addition. The base deprotonates the ketone at the only available alpha-carbon (C6).
    # The resulting enolate attacks the beta-carbon of ethyl acrylate.
    # When naming the product, the ester functional group has higher priority than the ketone.
    # The parent chain is ethyl propanoate. The cyclohexyl ring is a substituent at position 3.
    # Numbering the substituent ring from the point of attachment (original C6) gives:
    # - methyl at C1 (new)
    # - oxo at C2 (new)
    # - ethyl and methyl at C3 (new)
    # This leads to the substituent name (3-ethyl-1,3-dimethyl-2-oxocyclohexyl).
    correct_product_A = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"

    # Reaction B: 1-nitropropane + (KOH, (E)-but-2-enenitrile, H2O)
    # This is a nitro-Michael addition. The base deprotonates the carbon alpha to the nitro group.
    # The resulting carbanion attacks the beta-carbon of the unsaturated nitrile.
    # The principal functional group is the nitrile (-CN), so it's C1.
    # The longest chain containing the nitrile is 6 carbons long (hexanenitrile).
    # Numbering from the nitrile gives a methyl group at C3 and a nitro group at C4.
    correct_product_B = "3-methyl-4-nitrohexanenitrile"

    # Step 2: Define the product names for each option in the question.
    options = {
        "A": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "B": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        },
        "C": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "D": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        }
    }

    # Step 3: Extract the selected answer choice from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect format: The final answer was not found in the required <<<X>>> format."

    selected_option_key = match.group(1)

    # Step 4: Retrieve the product names corresponding to the selected answer.
    selected_products = options.get(selected_option_key)
    if not selected_products:
        return f"Invalid option '{selected_option_key}' selected. Valid options are A, B, C, D."

    # Step 5: Compare the selected product names with the correct ones and generate error messages if needed.
    errors = []
    if selected_products["A"] != correct_product_A:
        error_reason = (
            f"Product A is incorrect. The selected answer claims A is '{selected_products['A']}', "
            f"but the correct product is '{correct_product_A}'. "
            "Reasoning: In Reaction A, deprotonation of the ketone can only occur at C6. The subsequent Michael addition and IUPAC naming "
            "(ester > ketone priority) results in ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate."
        )
        errors.append(error_reason)

    if selected_products["B"] != correct_product_B:
        error_reason = (
            f"Product B is incorrect. The selected answer claims B is '{selected_products['B']}', "
            f"but the correct product is '{correct_product_B}'. "
            "Reasoning: In Reaction B, the product is named as a hexanenitrile, as the longest carbon chain containing the principal nitrile group has 6 carbons. "
            "This leads to the name 3-methyl-4-nitrohexanenitrile."
        )
        errors.append(error_reason)

    # Step 6: Return the final verdict.
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# The user provided the following text as the answer to be checked.
llm_answer_to_check = """
Here is a step-by-step analysis of the two reactions to determine the major products.

### **Reaction A: 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK) ---> A**

1.  **Identify the Reaction Type:** This is a Michael addition reaction. The reactants are a ketone (which can form an enolate, the Michael donor), an α,β-unsaturated ester (the Michael acceptor), and a strong base (t-BuOK).

2.  **Formation of the Nucleophile (Enolate):** The base, potassium tert-butoxide (t-BuOK), will abstract an acidic α-proton from the ketone.
    *   The ketone, 2-ethyl-2,6-dimethylcyclohexan-1-one, has two α-carbons (carbons adjacent to the C=O group at C1): C2 and C6.
    *   **C2:** This carbon is substituted with an ethyl group and a methyl group. It is a quaternary carbon and has **no α-protons**.
    *   **C6:** This carbon is substituted with a methyl group and has **one α-proton**.
    *   Therefore, deprotonation can only occur at the C6 position. The base removes this proton to form the enolate, which is nucleophilic at C6.

3.  **Nucleophilic Attack:** The C6 of the enolate attacks the electrophilic β-carbon of the Michael acceptor, ethyl acrylate (CH₂=CH-COOEt). A new carbon-carbon bond is formed between C6 of the cyclohexanone ring and the β-carbon of the ethyl acrylate.

4.  **Structure and Naming of Product A:** The resulting product has an ethyl propanoate group attached to the C6 position of the original ketone. The options name the product as a derivative of ethyl propanoate, which is the correct IUPAC convention as the ester group has higher priority than the ketone.
    *   The parent chain is **ethyl propanoate**.
    *   The cyclohexanone ring is a substituent on C3 of the propanoate chain.
    *   To name the cyclohexyl substituent, we number the ring starting from the point of attachment. The original C6 of the ketone is now C1 of the substituent.
    *   The substituents on the ring are:
        *   A **methyl** group at the new C1 (from the original C6-methyl).
        *   An **oxo** group at the new C2 (from the original C1 carbonyl).
        *   An **ethyl** group and a **methyl** group at the new C3 (from the original C2 substituents).
    *   Combining these gives the substituent name: **(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)**.
    *   The full name for product A is **ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate**.

### **Reaction B: 1-nitropropane + (KOH, (E)-but-2-enenitrile, H₂O) ---> B**

1.  **Identify the Reaction Type:** This is also a Michael addition. The nitroalkane has acidic α-protons and acts as the Michael donor.

2.  **Formation of the Nucleophile (Nitronate):** The base, KOH, deprotonates the carbon α to the electron-withdrawing nitro group in 1-nitropropane (CH₃CH₂CH₂NO₂). This forms a resonance-stabilized carbanion (a nitronate), which is nucleophilic at the α-carbon.

3.  **Nucleophilic Attack:** The nucleophilic α-carbon of the nitronate attacks the electrophilic β-carbon of the Michael acceptor, (E)-but-2-enenitrile (CH₃-CH=CH-CN).

4.  **Structure and Naming of Product B:** A new C-C bond forms between the α-carbon of the 1-nitropropane and the β-carbon of the but-2-enenitrile.
    *   The resulting structure is: `CH₃-CH₂-CH(NO₂)-CH(CH₃)-CH₂-CN`.
    *   To name this molecule, we identify the principal functional group, which is the nitrile (-CN). The carbon of the nitrile group is C1.
    *   The longest carbon chain containing C1 has 6 carbons: `CN(1)-CH₂(2)-CH(CH₃)(3)-CH(NO₂)(4)-CH₂(5)-CH₃(6)`.
    *   The parent name is **hexanenitrile**.
    *   The substituents are a **methyl** group at C3 and a **nitro** group at C4.
    *   The full name for product B is **3-methyl-4-nitrohexanenitrile**.

### **Conclusion**

*   **Product A:** ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate
*   **Product B:** 3-methyl-4-nitrohexanenitrile

Comparing these correct names with the given options:
*   A) A = ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate, B = 2,3-dimethyl-4-nitrobutanenitrile (Incorrect A, Incorrect B)
*   B) A = ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate, B = 3-methyl-4-nitrohexanenitrile (**Correct A, Correct B**)
*   C) A = ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate, B = 2,3-dimethyl-4-nitrobutanenitrile (Correct A, Incorrect B)
*   D) A = ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate, B = 3-methyl-4-nitrohexanenitrile (Incorrect A, Correct B)

The only option that correctly identifies both major products is B.

<<<B>>>
"""

# print(check_answer_correctness(llm_answer_to_check))