import re

def check_cycloaddition_product():
    """
    Checks the correctness of the final answer for the given [4+2] cycloaddition question.
    The function verifies the product's structural name and its stereochemistry.
    """
    # --- Problem Definition & Provided Answer ---
    question = {
        "diene": "2,5-dimethylthiophene",
        "dienophile": "Furan-2,5-dione",
        "product_type": "EXO"
    }
    
    final_answer_from_llm = "C"
    
    options = {
        "A": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "C": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "D": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione"
    }

    # --- Step 1: Verify the selected option exists ---
    selected_name = options.get(final_answer_from_llm)
    if not selected_name:
        return f"Incorrect. The provided answer '{final_answer_from_llm}' is not one of the valid options (A, B, C, D)."

    # --- Step 2: Verify the product's core structure and connectivity ---
    # The diene (thiophene derivative) provides a sulfur bridge, named "epithio".
    # The dienophile (maleic anhydride) provides the "isobenzofuran-1,3-dione" core.
    
    # Check for the correct bridge type
    if "epoxy" in selected_name:
        return (f"Incorrect. The answer '{final_answer_from_llm}' is wrong because it describes a product with an 'epoxy' (oxygen) bridge. "
                f"The diene is 2,5-dimethylthiophene, so the product must have an 'epithio' (sulfur) bridge.")
    
    # Check for the correct base name
    if "isobenzofuran-1,3-dione" not in selected_name:
        return (f"Incorrect. The answer '{final_answer_from_llm}' is wrong because it has an incorrect base name. "
                f"The structure derived from the maleic anhydride dienophile should be named as an 'isobenzofuran-1,3-dione' derivative.")

    # --- Step 3: Verify the product's stereochemistry ---
    # The question asks for the EXO product. We need to check if the selected answer's
    # stereochemical descriptor matches the known configuration for the EXO isomer.
    
    known_exo_config = "(3aR,4S,7R,7aS)"
    known_endo_config = "(3aR,4R,7S,7aS)" # This corresponds to option A

    # Extract the stereochemical descriptor from the name
    stereo_match = re.match(r"^\(.*\)", selected_name)
    if not stereo_match:
        return f"Incorrect. Could not extract stereochemical descriptors from the name of option '{final_answer_from_llm}'."
    
    stereo_descriptor = stereo_match.group(0)

    if stereo_descriptor == known_exo_config:
        # The selected answer has the correct stereochemistry for the EXO product.
        # Since it also passed the structural checks, it is correct.
        return "Correct"
    elif stereo_descriptor == known_endo_config:
        return (f"Incorrect. The answer '{final_answer_from_llm}' is wrong because it has the stereochemistry '{stereo_descriptor}', "
                f"which corresponds to the ENDO product. The question specifically asks for the EXO product.")
    else:
        return (f"Incorrect. The stereochemistry of the answer '{final_answer_from_llm}' ('{stereo_descriptor}') does not match "
                f"the expected configuration for either the EXO or ENDO product.")

# Run the check
result = check_cycloaddition_product()
print(result)