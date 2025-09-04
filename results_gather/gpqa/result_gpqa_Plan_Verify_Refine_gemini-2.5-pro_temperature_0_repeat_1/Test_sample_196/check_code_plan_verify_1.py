import re

def check_correctness():
    """
    This function checks the correctness of the LLM's answer to a chemistry problem.
    It follows these steps:
    1.  Deduces the structure of the starting material (Compound X) from the provided spectral data.
    2.  Identifies the chemical reaction (reduction with Red P/HI) and determines the final product.
    3.  Compares this derived correct product with the LLM's chosen answer and reasoning.
    4.  Returns "Correct" if the LLM's answer and reasoning are sound, otherwise returns an explanation of the error.
    """

    # --- Data from the problem ---
    llm_answer_text = """
    My plan is to first determine the structure of the starting material, Compound X, by interpreting its IR and 1H NMR spectra. Then, I will identify the reaction that occurs with red phosphorus and HI and apply it to Compound X to find the final product.

    First, I will analyze the provided spectral data to deduce the structure of Compound X.

    *   **IR Analysis:**
        *   The very broad peak from `3400–2500 cm-1` is characteristic of the O-H stretch in a carboxylic acid, broadened by hydrogen bonding.
        *   The strong peak at `1720 cm-1` is characteristic of the C=O (carbonyl) stretch of a carboxylic acid.
        *   The peaks at `1610 cm-1` and `1450 cm-1` are typical for C=C stretches in an aromatic ring.
        *   **Conclusion from IR:** Compound X is an aromatic carboxylic acid.

    *   **1H NMR Analysis:**
        *   `10.5 ppm (bs, 1H)`: This broad singlet is definitive for the acidic proton of a carboxylic acid (-COOH).
        *   `8.0 ppm (d, 2H)` and `7.2 ppm (d, 2H)`: This pattern of two doublets, each integrating to 2H, is classic for a 1,4-disubstituted (para) benzene ring. The protons at 8.0 ppm are ortho to the electron-withdrawing -COOH group.
        *   This confirms we have a benzoic acid with a substituent at the para (4) position. Now let's identify the alkyl substituent from the remaining signals:
            *   `0.9 ppm (t, 3H)`: A triplet for 3H indicates a -CH₃ group next to a -CH₂ group.
            *   `1.4 ppm (d, 3H)`: A doublet for 3H indicates a -CH₃ group next to a -CH group.
            *   `1.7 ppm (m, 2H)`: A multiplet for 2H is a -CH₂- group, likely between the -CH₃ (at 0.9 ppm) and the -CH group (at 2.9 ppm).
            *   `2.9 ppm (m, 1H)`: A multiplet for 1H is a -CH group attached to the benzene ring (benzylic position), a -CH₃ group, and a -CH₂- group.
        *   Assembling these fragments (`-CH-`, `-CH₃`, and `-CH₂CH₃`) gives a -CH(CH₃)(CH₂CH₃) group, which is a **sec-butyl** group.

    *   **Structure of Compound X:** Based on the combined spectral data, Compound X is **4-(sec-butyl)benzoic acid**.

    Now, I will determine the effect of the reaction.
    *   **Reaction:** The reaction of a carboxylic acid with red phosphorus and hydroiodic acid (HI) is a powerful reduction. This reaction reduces the carboxylic acid group (-COOH) completely to a methyl group (-CH₃).
    *   **Transformation:**
        *   Starting Material: 4-(sec-butyl)benzoic acid
        *   The -COOH group on the ring is reduced to a -CH₃ group.
        *   The sec-butyl group and the aromatic ring are unaffected.
    *   **Final Product:** The final product is **1-(sec-butyl)-4-methylbenzene**.

    Comparing this result with the given options:
    A) 4-(sec-butyl)benzoic acid - This is the starting material.
    B) 2-(4-ethylphenyl)propanoic acid - Incorrect structure.
    C) 1-isobutyl-4-methylbenzene - This has an isobutyl group, not a sec-butyl group.
    D) 1-(sec-butyl)-4-methylbenzene - This matches our deduced product.

    <<<D>>>
    """
    
    options = {
        "A": "4-(sec-butyl)benzoic acid",
        "B": "2-(4-ethylphenyl)propanoic acid",
        "C": "1-isobutyl-4-methylbenzene",
        "D": "1-(sec-butyl)-4-methylbenzene"
    }

    # --- Step 1: Deduce the correct answer from the problem statement ---
    # IR (3400-2500, 1720) + NMR (10.5 ppm) -> Aromatic Carboxylic Acid
    # NMR (8.0, 7.2 ppm, 2H doublets) -> 1,4-disubstituted benzene
    # NMR (0.9, 1.4, 1.7, 2.9 ppm) -> sec-butyl group: -CH(CH3)(CH2CH3)
    # So, Compound X is 4-(sec-butyl)benzoic acid.
    correct_starting_material = "4-(sec-butyl)benzoic acid"

    # Reaction: Red P/HI reduces -COOH to -CH3.
    # Product: 1-(sec-butyl)-4-methylbenzene.
    correct_final_product = "1-(sec-butyl)-4-methylbenzene"
    correct_option_letter = "D"

    # --- Step 2: Check the LLM's answer and reasoning ---
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "The answer is malformed. It does not contain a final choice in the format <<<X>>>."
    
    llm_choice_letter = match.group(1)
    
    # Check if the LLM's reasoning for the starting material is correct
    if correct_starting_material not in llm_answer_text:
        return f"The reasoning is flawed. The starting material, Compound X, was not correctly identified as '{correct_starting_material}'."

    # Check if the LLM's reasoning for the reaction is correct
    if "reduces the carboxylic acid group (-COOH) completely to a methyl group (-CH₃)" not in llm_answer_text:
        return "The reasoning is flawed. The function of the reagents (Red P/HI) to reduce a carboxylic acid to a methyl group was not correctly identified."

    # Check if the LLM's final product identification is correct
    if correct_final_product not in llm_answer_text:
        return f"The reasoning is flawed. The final product was not correctly identified as '{correct_final_product}'."

    # Check if the final chosen option is correct
    if llm_choice_letter == correct_option_letter:
        return "Correct"
    else:
        return f"The final answer is incorrect. The correct option is {correct_option_letter} ('{correct_final_product}'), but the answer provided was {llm_choice_letter} ('{options[llm_choice_letter]}')."

# To display the result of the check, we call the function.
# In a real environment, this function would be executed.
# For this demonstration, we print its return value.
print(check_correctness())