import sys
import io

def check_diels_alder_product():
    """
    This function checks the correctness of the proposed answer for a Diels-Alder reaction question.
    It verifies two main constraints:
    1. The basic chemical structure (connectivity) based on the reactants.
    2. The specific stereochemistry (EXO vs. ENDO) requested by the question.
    """
    # --- Define the problem constraints and options ---
    question = {
        "diene": "2,5-dimethylthiophene",
        "dienophile": "Furan-2,5-dione",
        "product_type": "EXO"
    }

    options = {
        'A': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'B': "(3aR,4S,7R,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione",
        'C': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epoxybenzo[c]thiophene-1,3-dione",
        'D': "(3aR,4R,7S,7aS)-4,7-dimethyl-3a,4,7,7a-tetrahydro-4,7-epithioisobenzofuran-1,3-dione"
    }

    # The final answer provided by the LLM to be checked.
    proposed_answer_key = 'B'
    
    # --- Verification Logic ---

    # Constraint 1: Check the product's core structure based on reactants.
    # The diene is thiophene-based, so the bridge must be sulfur ("epithio").
    # An oxygen bridge ("epoxy") would come from a furan diene.
    # The dienophile is maleic anhydride, which forms an "isobenzofuran-1,3-dione" base structure.
    required_structural_keywords = ["epithio", "isobenzofuran-1,3-dione"]
    
    # Filter options that have the correct core structure.
    structurally_correct_options = []
    for key, name in options.items():
        if all(keyword in name for keyword in required_structural_keywords):
            structurally_correct_options.append(key)

    if proposed_answer_key not in structurally_correct_options:
        if "epoxy" in options[proposed_answer_key]:
            return f"Incorrect. The proposed answer {proposed_answer_key} describes a product with an 'epoxy' (oxygen) bridge. The diene is 2,5-dimethylthiophene, which must form an 'epithio' (sulfur) bridge."
        else:
            return f"Incorrect. The proposed answer {proposed_answer_key} has an incorrect base name. The product of a Diels-Alder reaction with maleic anhydride should be named as an 'isobenzofuran-1,3-dione' derivative."

    # Constraint 2: Check the stereochemistry.
    # The question asks for the EXO product.
    # For this class of reactions, the stereochemical descriptors are well-established.
    # EXO adduct: (3aR,4S,7R,7aS) or its enantiomer (3aS,4R,7S,7aR)
    # ENDO adduct: (3aR,4R,7S,7aS) or its enantiomer (3aS,4S,7R,7aR)
    
    exo_stereochem_descriptor = "(3aR,4S,7R,7aS)"
    endo_stereochem_descriptor = "(3aR,4R,7S,7aS)"

    if question["product_type"] != "EXO":
        return "Checker Error: The checker is configured to check for an EXO product, but the question asks for something else."

    proposed_answer_name = options[proposed_answer_key]
    
    if proposed_answer_name.startswith(exo_stereochem_descriptor):
        # The proposed answer has the correct structure and the correct EXO stereochemistry.
        # Let's ensure no other option is also correct.
        # Option D has the correct structure but ENDO stereochemistry.
        if options['D'].startswith(endo_stereochem_descriptor):
             return "Correct"
        else:
             return "Checker Error: The stereochemistry of option D is not as expected."
    elif proposed_answer_name.startswith(endo_stereochem_descriptor):
        return f"Incorrect. The proposed answer {proposed_answer_key} has the stereochemistry '{endo_stereochem_descriptor}', which corresponds to the ENDO product. The question asks for the EXO product."
    else:
        return f"Incorrect. The stereochemistry of the proposed answer {proposed_answer_key} does not match the known patterns for EXO or ENDO products."

# Execute the check and print the result
result = check_diels_alder_product()
print(result)