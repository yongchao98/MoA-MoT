def check_diels_alder_product_correctness():
    """
    This function programmatically checks the correctness of the provided answer for the
    Diels-Alder reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The check is performed by applying established chemical principles of stereoselectivity
    and verifying if the chosen answer aligns with the predicted major product.
    """

    # --- Problem Definition & LLM's Answer ---
    llm_provided_answer = "C"
    options = {
        "A": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # --- Verification using Chemical Principles ---

    # Principle 1: Endo/Exo Selectivity (The Endo Rule)
    # In Diels-Alder reactions involving cyclic dienes like cyclopentadiene, the 'endo'
    # adduct is the kinetically favored major product due to stabilizing secondary orbital interactions.
    # The 'exo' product is the minor product.
    # Based on IUPAC nomenclature for this bicyclic system:
    # - The stereochemistry (3aR,4R,7S,7aS) corresponds to the 'endo' adduct.
    # - The stereochemistry (3aR,4S,7R,7aS) corresponds to the 'exo' adduct.
    
    endo_stereochem_pattern = "(3aR,4R,7S,7aS"
    endo_candidates = {key for key, name in options.items() if endo_stereochem_pattern in name}
    
    # The major product must be an 'endo' adduct. We filter the options accordingly.
    if endo_candidates != {'A', 'C'}:
        return (f"Constraint check failed: The 'endo rule' dictates that the major product must have 'endo' "
                f"stereochemistry. This corresponds to options A and C. The logic is flawed if other "
                f"options are considered at this stage.")

    # Principle 2: Syn/Anti Facial Selectivity
    # The dienophile can attack from the same face as the C5-substituent ('syn' attack)
    # or from the opposite face ('anti' attack).
    # The substituent is Fluorine (F), a small, highly electronegative, electron-withdrawing group (EWG).
    # For a C5-EWG on cyclopentadiene, 'syn' attack is electronically favored over 'anti' attack, which would be favored by sterics.
    # Electronic effects typically dominate for fluorine, so the major product results from 'syn' attack.
    
    # The stereochemistry at the C8 position (bearing the fluorine) distinguishes syn/anti products.
    # For this specific bicyclic system, the relative descriptors are:
    # - 'syn' relationship (fluorine and anhydride on the same side) corresponds to the (8r) descriptor.
    # - 'anti' relationship (fluorine and anhydride on opposite sides) corresponds to the (8s) descriptor.
    
    syn_descriptor_pattern = ",8r)"
    
    # From the 'endo' candidates, we now select the one resulting from 'syn' attack.
    major_product_key = None
    for key in endo_candidates:
        if syn_descriptor_pattern in options[key]:
            major_product_key = key
            break  # Found the syn-endo product

    # --- Final Conclusion ---
    
    # Check if a valid major product was identified based on the principles.
    if major_product_key is None:
        return ("Constraint check failed: Could not identify a 'syn-endo' product among the options. "
                "The expected product should have 'endo' stereochemistry and a 'syn' (8r) fluorine.")

    # Check if the LLM's answer matches the logically derived major product.
    if llm_provided_answer == major_product_key:
        return "Correct"
    else:
        return (f"The answer is incorrect. "
                f"Based on established chemical principles, the major product is the 'syn-endo' adduct. "
                f"Step 1 (Endo Rule) narrows the choices to {endo_candidates}. "
                f"Step 2 (Syn-Attack Rule) selects the product with the 'syn' (8r) stereocenter from these candidates. "
                f"This logic correctly identifies option '{major_product_key}' as the major product. "
                f"The provided answer was '{llm_provided_answer}'.")

# To run the check, you would execute the function:
result = check_diels_alder_product_correctness()
print(result)