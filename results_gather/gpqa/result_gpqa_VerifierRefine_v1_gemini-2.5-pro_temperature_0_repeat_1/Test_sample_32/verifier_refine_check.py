import re

def check_diels_alder_product():
    """
    Checks the correctness of the selected answer for the Diels-Alder reaction
    between 2,5-dimethylthiophene and maleic anhydride.
    """
    # --- Problem Definition ---
    question = "Identify the EXO product of the following [4+2] cycloaddition reaction. 2,5-dimethylthiophene + Furan-2,5-dione + Heat ---> ?"
    options = {
        "A": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "B": "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        "D": "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }
    llm_answer_choice = "B"
    
    # --- Verification ---
    
    selected_answer_name = options.get(llm_answer_choice)
    if not selected_answer_name:
        return f"Error: The provided answer choice '{llm_answer_choice}' is not a valid option."

    # 1. Check the bridge atom (Connectivity)
    # The diene is thiophene, so the bridge atom must be sulfur (S).
    # The correct prefix for a sulfur bridge is "epithio". "Epoxy" implies an oxygen bridge.
    if 'epoxy' in selected_answer_name:
        return (f"Incorrect. The answer '{llm_answer_choice}' is wrong because it describes a product with an 'epoxy' (oxygen) bridge. "
                "The diene is 2,5-dimethylthiophene, so the bridge atom must be sulfur, correctly named as 'epithio'.")

    if 'epithio' not in selected_answer_name:
        return (f"Incorrect. The answer '{llm_answer_choice}' is wrong because the product name does not contain the required 'epithio' "
                "prefix to describe the sulfur bridge formed from the thiophene diene.")

    # 2. Check the stereochemistry (EXO vs. ENDO)
    # The question specifically asks for the EXO product. We need to check if the
    # stereochemical descriptors match the known EXO adduct.
    
    # Define the stereodescriptors for the two enantiomers of the EXO product.
    exo_stereochem_patterns = [
        r'\(3aR,4S,7R,7aS\)',  # One enantiomer
        r'\(3aS,4R,7S,7aR\)'   # The other enantiomer (mirror image)
    ]

    # Define the stereodescriptors for the ENDO product for comparison.
    endo_stereochem_patterns = [
        r'\(3aR,4R,7S,7aS\)',  # One enantiomer (matches option D)
        r'\(3aS,4S,7R,7aR\)'   # The other enantiomer
    ]

    is_exo = any(re.search(pattern, selected_answer_name) for pattern in exo_stereochem_patterns)
    is_endo = any(re.search(pattern, selected_answer_name) for pattern in endo_stereochem_patterns)

    if is_endo:
        return (f"Incorrect. The answer '{llm_answer_choice}' is wrong because it describes the ENDO product. "
                "The question specifically asks for the EXO product, which is a different stereoisomer.")

    if not is_exo:
        return (f"Incorrect. The answer '{llm_answer_choice}' is wrong because its stereochemical descriptor "
                "does not match either of the two possible enantiomers of the EXO product.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_diels_alder_product()
print(result)
