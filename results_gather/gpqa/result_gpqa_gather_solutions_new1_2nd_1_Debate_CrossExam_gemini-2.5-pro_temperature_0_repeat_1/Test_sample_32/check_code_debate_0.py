import re

def check_correctness():
    """
    Checks the correctness of the provided LLM answer for the chemical reaction question.

    The function verifies the answer against three key constraints:
    1.  **Bridge Type**: The diene is 2,5-dimethylthiophene, so the product must have a sulfur bridge, named "epithio".
    2.  **Base Name**: The dienophile is maleic anhydride, which forms an "isobenzofuran-1,3-dione" core in the product.
    3.  **Stereochemistry**: The question asks for the EXO product. Based on Cahn-Ingold-Prelog rules, the EXO isomer has the stereochemical descriptor (3aR,4S,7R,7aS).
    """
    # Options as defined in the original question prompt.
    options = {
        "A": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    selected_name = options.get(llm_answer)
    if not selected_name:
        return f"Invalid option '{llm_answer}' selected. Not one of A, B, C, D."

    # --- Constraint 1: Bridge Type ---
    # The diene (2,5-dimethylthiophene) contains sulfur, so the bridge must be "epithio".
    if "epoxy" in selected_name:
        return (f"Incorrect: The selected answer '{llm_answer}' describes a product with an 'epoxy' (oxygen) bridge. "
                "This is wrong because the diene was thiophene, which contains sulfur. The bridge must be 'epithio'.")
    if "epithio" not in selected_name:
        return (f"Incorrect: The selected answer '{llm_answer}' does not describe a product with the correct 'epithio' (sulfur) bridge, "
                "which is expected from the thiophene diene.")

    # --- Constraint 2: Base Name ---
    # The dienophile (maleic anhydride) forms an "isobenzofuran-1,3-dione" core.
    correct_base_name = "isobenzofuran-1,3-dione"
    if correct_base_name not in selected_name:
        return (f"Incorrect: The selected answer '{llm_answer}' has an incorrect base name. "
                f"The structure derived from maleic anhydride should be named '{correct_base_name}'.")

    # --- Constraint 3: Stereochemistry (EXO product) ---
    # The question asks for the EXO product. Chemical analysis shows the EXO isomer
    # has the stereodescriptor (3aR,4S,7R,7aS).
    correct_exo_stereochemistry = "(3aR,4S,7R,7aS)"
    endo_stereochemistry = "(3aR,4R,7S,7aS)"
    
    if not selected_name.startswith(correct_exo_stereochemistry):
        if selected_name.startswith(endo_stereochemistry):
            return (f"Incorrect: The selected answer '{llm_answer}' has the stereochemistry for the ENDO product. "
                    f"The question specifically asks for the EXO product, which has the stereochemistry '{correct_exo_stereochemistry}'.")
        else:
            return (f"Incorrect: The selected answer '{llm_answer}' does not have the correct stereochemistry for the EXO product. "
                    f"The correct stereodescriptor should be '{correct_exo_stereochemistry}'.")

    # If all constraints are satisfied
    return "Correct"

# Run the check and print the result
result = check_correctness()
print(result)