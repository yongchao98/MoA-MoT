import re

def check_answer():
    """
    Checks the correctness of the final answer for the Diels-Alder reaction question.
    
    The function verifies the answer against three main constraints:
    1.  The fundamental structure of the product (nomenclature).
    2.  The stereochemistry (EXO vs. ENDO).
    3.  The specific R/S configuration corresponding to the EXO product.
    """
    
    # The final consolidated answer to be checked
    final_answer_llm = "<<<A>>>"

    # --- Define the problem constraints and correct properties ---

    # Constraint 1: Reaction type and reactants
    # Diene: 2,5-dimethylthiophene -> provides a sulfur bridge ("epithio")
    # Dienophile: Furan-2,5-dione -> forms an "isobenzofuran-1,3-dione" core
    correct_base_name = "epithioisobenzofuran-1,3-dione"

    # Constraint 2: Desired product stereoisomer
    # The question asks for the "EXO" product.
    # The "Heat" condition favors the thermodynamically more stable product, which is EXO.
    
    # Constraint 3: Correct R/S configuration for the EXO product
    # Based on Cahn-Ingold-Prelog rules and analogy to similar reactions (e.g., furan + maleic anhydride),
    # the EXO adduct has the (3aR,4S,7R,7aS) configuration.
    # The ENDO adduct has the (3aR,4R,7S,7aS) configuration.
    correct_exo_stereochem = "(3aR,4S,7R,7aS)"
    
    # --- Define the options as presented in the final consolidated answer ---
    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'C': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'D': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # --- Verification Logic ---

    # 1. Extract the letter from the LLM's final answer
    match = re.search(r'<<<([A-D])>>>', final_answer_llm)
    if not match:
        return f"Invalid answer format. The answer should be in the format '<<<X>>>' where X is A, B, C, or D. Received: {final_answer_llm}"
    
    chosen_letter = match.group(1)
    chosen_answer_text = options[chosen_letter]

    # 2. Check the product's core structure (base name)
    if correct_base_name not in chosen_answer_text:
        return (f"Incorrect product structure. The chosen answer '{chosen_letter}' describes a product with an incorrect core name. "
                f"The reaction between 2,5-dimethylthiophene and furan-2,5-dione should yield an '{correct_base_name}' derivative. "
                f"The chosen answer describes a '{chosen_answer_text.split('-')[-1]}'.")

    # 3. Check the stereochemistry for the EXO product
    if correct_exo_stereochem not in chosen_answer_text:
        # This means the chosen answer has the correct base name but wrong stereochemistry (it's the ENDO product)
        actual_stereochem = chosen_answer_text.split('-')[0]
        return (f"Incorrect stereochemistry. The chosen answer '{chosen_letter}' describes the ENDO product, not the EXO product. "
                f"The question asks for the EXO product, which has the '{correct_exo_stereochem}' configuration. "
                f"The chosen answer has the '{actual_stereochem}' configuration.")

    # 4. If all checks pass, the answer is correct
    return "Correct"

# Execute the check and print the result
result = check_answer()
print(result)