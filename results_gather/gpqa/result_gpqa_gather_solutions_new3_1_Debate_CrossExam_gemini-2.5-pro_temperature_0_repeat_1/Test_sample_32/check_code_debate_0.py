import re

def check_diels_alder_product():
    """
    Checks the correctness of the selected answer for the Diels-Alder reaction.
    """
    # --- Problem Definition ---
    question = {
        "diene": "2,5-dimethylthiophene",
        "dienophile": "Furan-2,5-dione",
        "product_type": "EXO"
    }
    
    options = {
        'A': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'C': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'D': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }
    
    # The final answer to be checked
    final_answer = 'C'
    
    # --- Verification Logic ---
    
    selected_name = options.get(final_answer)
    if not selected_name:
        return f"Error: The provided answer '{final_answer}' is not one of the options."

    # 1. Check Connectivity Constraint: The bridge must be "epithio" (from thiophene).
    if "epoxy" in selected_name:
        return (f"Incorrect. The answer '{final_answer}' is wrong because its name contains 'epoxy', "
                f"indicating an oxygen bridge. The diene was {question['diene']}, which contains sulfur, "
                f"so the product must have an 'epithio' bridge.")

    if "epithioisobenzofuran-1,3-dione" not in selected_name:
        return (f"Incorrect. The answer '{final_answer}' has an incorrect base name. The product of a "
                f"thiophene derivative and maleic anhydride should be named as a derivative of "
                f"'epithioisobenzofuran-1,3-dione'.")

    # 2. Check Stereochemistry Constraint: The product must be EXO.
    # Based on chemical literature and analysis, we can assign stereodescriptors to EXO/ENDO.
    known_exo_config = "3aR,4S,7R,7aS"
    known_endo_config = "3aR,4R,7S,7aS"
    
    # Extract the stereodescriptor from the selected name
    match = re.match(r'\((.*?)\)', selected_name)
    if not match:
        return f"Incorrect. The name for option '{final_answer}' is missing stereochemical information."
    
    actual_config = match.group(1)

    if question["product_type"] == "EXO":
        if actual_config == known_exo_config:
            # The configuration matches the known EXO configuration.
            pass
        elif actual_config == known_endo_config:
            return (f"Incorrect. The answer '{final_answer}' is wrong because it describes the ENDO product. "
                    f"Its stereochemistry '{actual_config}' corresponds to the ENDO isomer, but the question "
                    f"asked for the EXO product.")
        else:
            # This case handles potential typos or other incorrect stereodescriptors.
            return (f"Incorrect. The stereochemistry '{actual_config}' for answer '{final_answer}' does not "
                    f"match the known configuration for the EXO product.")
    else:
        # This part would handle a question asking for the ENDO product.
        return "Checker Error: This check is configured for an EXO product question."

    # If all constraints are satisfied
    return "Correct"

# Execute the check and print the result
result = check_diels_alder_product()
print(result)