import re

def check_chemistry_answer(final_answer_text: str) -> str:
    """
    Checks the correctness of the final answer for the given organic chemistry problem.

    The function simulates the reasoning process:
    1. Analyzes spectral data to identify the starting material (Compound X).
    2. Determines the product of the reaction with red phosphorus and HI.
    3. Matches the product to the given options and compares with the provided answer.

    Args:
        final_answer_text: The full text of the final answer, including the <<<A>>> part.

    Returns:
        A string indicating "Correct" or the reason for the incorrectness.
    """

    # --- Step 0: Define the problem's options and extract the given answer ---
    
    # The options as listed in the final analysis provided by the user.
    options = {
        'A': '1-(sec-butyl)-4-methylbenzene',
        'B': '1-isobutyl-4-methylbenzene',
        'C': '4-(sec-butyl)benzoic acid',
        'D': '2-(4-ethylphenyl)propanoic acid'
    }

    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "The final answer is not in the required format '<<<A>>>', '<<<B>>>', '<<<C>>>', or '<<<D>>>'."
    given_answer = match.group(1)

    # --- Step 1: Identify the starting material (Compound X) from spectral data ---
    
    # This section codifies the interpretation of the spectra.
    # IR Data: 3400–2500 (broad), 1720, 1610, 1450 cm-1
    # 1H NMR Data: 10.5(bs,1H), 8.0(d,2H), 7.2(d,2H), 2.9(m,1H), 1.7(m,2H), 1.4(d,3H), 0.9(t,3H)
    
    # Check for carboxylic acid features
    has_cooh_ir = True  # Based on 3400-2500 and 1720 cm-1 peaks
    has_cooh_nmr = True # Based on 10.5 ppm peak
    
    # Check for aromatic features
    has_aromatic_ir = True # Based on 1610, 1450 cm-1 peaks
    has_para_ring_nmr = True # Based on two doublets at 8.0 and 7.2 ppm, each 2H

    # Check for alkyl group features (sec-butyl)
    has_sec_butyl_nmr = True # Based on the combination of t, d, and two m signals with correct integrations

    if not (has_cooh_ir and has_cooh_nmr and has_aromatic_ir and has_para_ring_nmr and has_sec_butyl_nmr):
        return "Reason: The logical deduction of the starting material from the spectral data is flawed. The data points to 4-(sec-butyl)benzoic acid."

    starting_material = "4-(sec-butyl)benzoic acid"

    # --- Step 2: Determine the product of the reaction ---

    # Reaction: Starting material + red phosphorus and HI
    # This is a known strong reduction that converts a carboxylic acid (-COOH) to a methyl group (-CH3).
    # The sec-butyl group and aromatic ring are unaffected.
    
    if starting_material == "4-(sec-butyl)benzoic acid":
        correct_product = "1-(sec-butyl)-4-methylbenzene"
    else:
        # This case should not be reached if Step 1 is correct.
        return "Reason: The starting material was not correctly identified, so the product cannot be determined."

    # --- Step 3: Match the product to the options and check the answer ---

    # Find which option corresponds to the correct product
    correct_option_letter = None
    for letter, name in options.items():
        if name == correct_product:
            correct_option_letter = letter
            break
    
    if correct_option_letter is None:
        return f"Reason: The deduced correct product '{correct_product}' does not match any of the provided options (A, B, C, D)."

    # Check if the given answer matches the derived correct option
    if given_answer == correct_option_letter:
        return "Correct"
    else:
        # Provide a detailed reason for the error
        if options.get(given_answer) == starting_material:
            return f"Reason: The provided answer '{given_answer}' is incorrect. It identifies the starting material, not the final product. The correct option is '{correct_option_letter}'."
        else:
            return f"Reason: The provided answer '{given_answer}' is incorrect. The correct product is '{correct_product}', which corresponds to option '{correct_option_letter}'."

# The user's provided text to check
final_answer_to_check = """
Here is a step-by-step analysis to determine the final answer.

### **Step 1: Determine the structure of the starting material, Compound X.**

This requires a careful analysis of the provided spectroscopic data.

*   **IR Spectrum Analysis:**
    *   The very broad absorption from **3400–2500 cm⁻¹** is a classic signature for the O-H stretch of a hydrogen-bonded carboxylic acid.
    *   The strong absorption at **1720 cm⁻¹** is characteristic of the C=O (carbonyl) stretch of a carboxylic acid.
    *   The peaks at **1610 cm⁻¹** and **1450 cm⁻¹** indicate the presence of C=C bonds within an aromatic ring.
    *   **Conclusion:** The IR data strongly suggests that Compound X is an aromatic carboxylic acid.

*   **¹H NMR Spectrum Analysis:**
    *   **10.5 ppm (bs, 1H):** This broad singlet in the far downfield region confirms the presence of the acidic proton of a carboxylic acid (-COOH).
    *   **8.0 ppm (d, 2H) and 7.2 ppm (d, 2H):** This pattern of two doublets, each integrating to two protons, is a hallmark of a 1,4-disubstituted (para) benzene ring. The protons at 8.0 ppm are ortho to the electron-withdrawing -COOH group, and the protons at 7.2 ppm are ortho to the other substituent.
    *   **Alkyl Group Analysis:** The remaining signals must be assembled to identify the second substituent.
        *   `0.9 ppm (t, 3H)` -> A -CH₃ group next to a -CH₂ group.
        *   `1.4 ppm (d, 3H)` -> A -CH₃ group next to a -CH group.
        *   `1.7 ppm (m, 2H)` -> A -CH₂- group.
        *   `2.9 ppm (m, 1H)` -> A -CH- group attached to the benzene ring.
        *   Assembling these pieces: The benzylic -CH- group is attached to a methyl group (-CH₃) and an ethyl group (-CH₂CH₃). This forms a **sec-butyl** group.
    *   **Conclusion:** Compound X is **4-(sec-butyl)benzoic acid**.

### **Step 2: Analyze the chemical reaction.**

*   **Reactant:** 4-(sec-butyl)benzoic acid (Compound X).
*   **Reagents:** Red phosphorus (P) and hydroiodic acid (HI).
*   **Transformation:** This combination of reagents is a very powerful reducing agent. Its specific function here is the complete reduction of the carboxylic acid group (-COOH) to a methyl group (-CH₃). The aromatic ring and the sec-butyl alkyl group are stable under these conditions and do not react or rearrange.

### **Step 3: Identify the final product.**

*   The reaction starts with 4-(sec-butyl)benzoic acid.
*   The -COOH group is reduced to a -CH₃ group.
*   The final product is therefore **1-(sec-butyl)-4-methylbenzene**.

### **Step 4: Match the final product with the given options.**

The options provided in the question are:
A) 1-(sec-butyl)-4-methylbenzene
B) 1-isobutyl-4-methylbenzene
C) 4-(sec-butyl)benzoic acid
D) 2-(4-ethylphenyl)propanoic acid

*   Our deduced final product is **1-(sec-butyl)-4-methylbenzene**, which corresponds to **Option A**.
*   Option B is incorrect because the alkyl group is sec-butyl, not isobutyl.
*   Option C is incorrect because it is the starting material (Compound X), not the final product.
*   Option D is an incorrect structure.

The final answer is A.

<<<A>>>
"""

# Run the check
result = check_chemistry_answer(final_answer_to_check)
print(result)