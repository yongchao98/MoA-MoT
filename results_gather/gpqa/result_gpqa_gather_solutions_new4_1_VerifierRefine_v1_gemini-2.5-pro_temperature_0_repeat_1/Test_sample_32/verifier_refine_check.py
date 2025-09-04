import re

def check_cycloaddition_answer():
    """
    Checks the correctness of the provided answer for the Diels-Alder reaction.
    """
    # The final answer provided by the LLM to be checked.
    llm_final_answer = 'C'

    # Define the options from the question.
    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'B': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'C': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'D': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }

    # --- Step 1: Verify the product's core structure from the reactants ---
    # Diene: 2,5-dimethylthiophene -> Bridge must be sulfur ("epithio").
    # Dienophile: Furan-2,5-dione -> Base name is a derivative of "isobenzofuran-1,3-dione".
    
    expected_bridge_keyword = "epithio"
    incorrect_bridge_keyword = "epoxy"
    expected_base_name_keyword = "isobenzofuran-1,3-dione"

    selected_option_text = options.get(llm_final_answer)

    if not selected_option_text:
        return f"Invalid answer choice '{llm_final_answer}'. Please choose from A, B, C, or D."

    # Check for correct bridge type
    if incorrect_bridge_keyword in selected_option_text:
        return (f"Incorrect. The answer '{llm_final_answer}' is wrong because it describes an 'epoxy' (oxygen) bridge. "
                f"The diene is 2,5-dimethylthiophene, which must form an 'epithio' (sulfur) bridge.")

    # Check for correct base name
    if expected_base_name_keyword not in selected_option_text:
        return (f"Incorrect. The answer '{llm_final_answer}' is wrong because it uses an incorrect base name. "
                f"The product should be named as a derivative of '{expected_base_name_keyword}'.")

    # --- Step 2: Verify the stereochemistry for the EXO product ---
    # The question asks for the EXO product. Based on established chemical principles,
    # the stereochemistry for the EXO and ENDO adducts can be determined.
    
    # Stereochemistry for the EXO product (one enantiomer)
    exo_stereochemistry = "(3aR,4S,7R,7aS)"
    
    # Stereochemistry for the ENDO product (one enantiomer)
    endo_stereochemistry = "(3aR,4R,7S,7aS)"

    # Extract the stereochemistry from the selected answer's text.
    match = re.match(r'(\(.*\))', selected_option_text)
    if not match:
        return f"Error: Could not parse the stereochemical descriptor from the answer text for option {llm_final_answer}."
    
    answer_stereochemistry = match.group(1)

    # Check if the answer's stereochemistry matches the required EXO product.
    if answer_stereochemistry == exo_stereochemistry:
        # This would be the correct case.
        return "Correct"
    elif answer_stereochemistry == endo_stereochemistry:
        # This is the case for the provided answer 'C'.
        return (f"Incorrect. The answer '{llm_final_answer}' identifies the product with stereochemistry {answer_stereochemistry}. "
                f"This configuration corresponds to the ENDO adduct, which is the kinetic product. "
                f"The question asks for the EXO product (the thermodynamic product, favored by heat), "
                f"which has the stereochemistry {exo_stereochemistry} (as seen in option D).")
    else:
        return (f"Incorrect. The stereochemistry {answer_stereochemistry} for option {llm_final_answer} does not match "
                f"the known configuration for either the EXO or ENDO product.")

# Run the check
result = check_cycloaddition_answer()
print(result)