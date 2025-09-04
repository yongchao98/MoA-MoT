import re

def check_cycloaddition_answer():
    """
    This function checks the correctness of the final answer for the given
    [4+2] cycloaddition reaction question by verifying two key constraints:
    1. The presence of the correct bridge atom (sulfur -> "epithio").
    2. The correct stereochemistry for the requested EXO product.
    """
    
    # --- Problem Definition ---
    # The final answer provided by the LLM to be checked.
    final_answer_key = 'A'

    # A dictionary mapping the option keys to their full IUPAC names.
    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'C': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'D': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # --- Constraint Definitions ---
    # The reactant is thiophene (sulfur), so the bridge must be "epithio".
    required_bridge_type = "epithio"
    
    # The question asks for the EXO product. For this molecule, the EXO isomer
    # has the stereochemical descriptor (3aR,4S,7R,7aS).
    required_stereochemistry_for_exo = "(3aR,4S,7R,7aS)"

    # --- Verification ---
    chosen_product_name = options.get(final_answer_key)

    if not chosen_product_name:
        return f"Error: The provided answer key '{final_answer_key}' is not a valid option."

    # 1. Check the structural constraint (bridge type).
    if required_bridge_type not in chosen_product_name:
        return (f"Incorrect. The chosen answer '{final_answer_key}' violates the structural constraint. "
                f"The reactant is thiophene (contains sulfur), so the product must have a sulfur bridge, "
                f"named with the prefix '{required_bridge_type}'. The chosen answer describes an 'epoxy' (oxygen) bridge.")

    # 2. Check the stereochemical constraint (EXO vs. ENDO).
    # We extract the stereochemistry part from the beginning of the name.
    match = re.match(r'\((.*?)\)', chosen_product_name)
    if not match:
        return f"Incorrect. The chosen answer '{final_answer_key}' is missing the stereochemical descriptor."
    
    chosen_stereochemistry = f"({match.group(1)})"

    if chosen_stereochemistry != required_stereochemistry_for_exo:
        return (f"Incorrect. The chosen answer '{final_answer_key}' violates the stereochemical constraint. "
                f"The question asks for the EXO product. The correct stereochemical descriptor for the EXO isomer is "
                f"'{required_stereochemistry_for_exo}'. The chosen answer has the descriptor '{chosen_stereochemistry}', "
                f"which corresponds to the ENDO isomer.")

    # If both constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_cycloaddition_answer()
print(result)