import re

def check_chemistry_answer():
    """
    Checks the correctness of the final answer for the given organic chemistry problem.

    The function verifies the products of two Michael addition reactions by:
    1. Defining the correct IUPAC names for the products based on established chemical principles.
    2. Defining the multiple-choice options provided in the question.
    3. Parsing the final answer given by the LLM.
    4. Comparing the products in the chosen option with the correct product names.
    5. Returning "Correct" if they match, or a detailed reason for the mismatch.
    """
    # The final answer text provided by the LLM to be checked.
    llm_answer_text = """
    To decide the final answer, a careful, step-by-step analysis of each reaction is required, followed by a precise matching of the resulting product names to the given options. Many of the provided candidate answers contain correct chemical reasoning but fail in the final step of matching their results to the options, leading to incorrect conclusions.

    ### **Analysis of Reaction A**

    **Reaction:** 2-ethyl-2,6-dimethylcyclohexan-1-one + ethyl acrylate (t-BuOK) ---> A

    1.  **Reaction Type:** This is a Michael addition. The ketone is the Michael donor, and ethyl acrylate is the Michael acceptor. The strong, bulky base (t-BuOK) generates the nucleophilic enolate from the ketone.

    2.  **Careful Point 1: Regioselectivity of Enolate Formation.** The most critical step is determining where the ketone is deprotonated. The ketone has two α-carbons (carbons adjacent to the carbonyl): C2 and C6.
        *   **C2** is a quaternary carbon (bonded to an ethyl and a methyl group). It has **no α-protons** and cannot be deprotonated.
        *   **C6** is a tertiary carbon (bonded to a methyl group and one hydrogen). It has **one α-proton**.
        *   **Conclusion:** The reaction is perfectly regioselective. The base can only abstract the proton from C6, forming a single, specific enolate.

    3.  **Product Structure:** The nucleophilic C6 of the enolate attacks the β-carbon of ethyl acrylate. A new C-C bond forms, attaching a `-CH₂CH₂COOEt` group to the C6 position of the original ketone.

    4.  **Careful Point 2: IUPAC Naming of Product A.** The ester functional group has higher priority than the ketone, so the product is named as a substituted ethyl propanoate.
        *   **Parent Chain:** `ethyl propanoate`.
        *   **Substituent:** The entire cyclohexanone ring is a substituent on the C3 position of the propanoate chain.
        *   **Naming the Substituent:** The ring is numbered starting from the point of attachment (the original C6) as the new C1. The numbering proceeds to give the ketone's oxo group the lowest possible number, which is C2.
            *   **New C1 (old C6):** Has a **methyl** group.
            *   **New C2 (old C1):** Has the **oxo** group.
            *   **New C3 (old C2):** Has an **ethyl** group and a **methyl** group.
        *   The correct substituent name is **(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)**.
        *   The full name for Product A is **ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate**.

    5.  **Constraint 1:** This name for Product A matches the name given in options **C** and **D**. Therefore, options A and B are eliminated.

    ### **Analysis of Reaction B**

    **Reaction:** 1-nitropropane + (KOH, (E)-but-2-enenitrile, H₂O) ---> B

    1.  **Reaction Type:** This is also a Michael addition. The nitroalkane is the Michael donor.

    2.  **Nucleophile Formation:** The base (KOH) deprotonates the carbon α to the electron-withdrawing nitro group in 1-nitropropane, forming the nucleophile `[CH₃CH₂CH-NO₂]⁻`.

    3.  **Product Structure:** The nucleophile attacks the β-carbon of (E)-but-2-enenitrile. The resulting structure is `CH₃CH₂-CH(NO₂)-CH(CH₃)-CH₂-CN`.

    4.  **Careful Point 3: IUPAC Naming of Product B.** The nitrile (-CN) is the principal functional group, and its carbon is C1.
        *   **Parent Chain:** The longest carbon chain containing the nitrile carbon has **six carbons**. The parent name is **hexanenitrile**. (This eliminates names based on butanenitrile).
        *   **Numbering and Substituents:** Numbering the chain from the nitrile carbon: `CN(1)-CH₂(2)-CH(CH₃)(3)-CH(NO₂)(4)-CH₂(5)-CH₃(6)`.
        *   The substituents are a **methyl** group at C3 and a **nitro** group at C4.
        *   The full name for Product B is **3-methyl-4-nitrohexanenitrile**.

    5.  **Constraint 2:** This name for Product B matches the name given in options **B** and **D**. Therefore, options A and C are eliminated.

    ### **Final Conclusion**

    *   **Constraint from Reaction A:** The correct option must be **C** or **D**.
    *   **Constraint from Reaction B:** The correct option must be **B** or **D**.

    The only option that satisfies both constraints is **D**.

    <<<D>>>
    """

    # Step 1: Define the correct product names based on chemical principles.
    correct_product_A = "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate"
    correct_product_B = "3-methyl-4-nitrohexanenitrile"

    # Step 2: Define the options from the question.
    options = {
        "A": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "B": {
            "A": "ethyl 3-(3-ethyl-3,5-dimethyl-4-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        },
        "C": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "2,3-dimethyl-4-nitrobutanenitrile"
        },
        "D": {
            "A": "ethyl 3-(3-ethyl-1,3-dimethyl-2-oxocyclohexyl)propanoate",
            "B": "3-methyl-4-nitrohexanenitrile"
        }
    }

    # Step 3: Extract the final answer from the provided text.
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Failure: Could not find the final answer in the format <<<A>>>, <<<B>>>, etc."
    
    chosen_option_key = match.group(1)
    chosen_option_products = options[chosen_option_key]

    # Step 4: Compare the chosen option's products with the correct ones.
    product_A_in_choice = chosen_option_products["A"]
    product_B_in_choice = chosen_option_products["B"]

    is_A_correct = (product_A_in_choice == correct_product_A)
    is_B_correct = (product_B_in_choice == correct_product_B)

    # Step 5 & 6: Return the result.
    if is_A_correct and is_B_correct:
        return "Correct"
    else:
        error_messages = [f"The final answer '{chosen_option_key}' is incorrect."]
        if not is_A_correct:
            error_messages.append(
                f"Reason for Product A: The name provided is '{product_A_in_choice}', but the correct name is '{correct_product_A}'. "
                "The error stems from the regioselectivity and IUPAC naming. The enolate can only form at the C6 position of the ketone (as C2 has no alpha-protons). "
                "When naming the resulting substituent, the point of attachment (the original C6) becomes C1, and the carbonyl group is at C2, not C4."
            )
        
        if not is_B_correct:
            error_messages.append(
                f"Reason for Product B: The name provided is '{product_B_in_choice}', but the correct name is '{correct_product_B}'. "
                "The error is in determining the parent chain length. The reaction between 1-nitropropane (3 carbons) and but-2-enenitrile (4 carbons) results in a product whose main chain, including the nitrile carbon, is 6 carbons long (hexanenitrile), not 4 carbons long (butanenitrile)."
            )
        
        return "\n".join(error_messages)

# Execute the check and print the result.
result = check_chemistry_answer()
print(result)