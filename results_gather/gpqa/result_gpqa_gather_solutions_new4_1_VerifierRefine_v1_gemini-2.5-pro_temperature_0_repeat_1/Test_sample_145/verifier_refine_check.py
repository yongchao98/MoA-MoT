import re

def check_diels_alder_product():
    """
    Checks the correctness of the predicted major product for the Diels-Alder reaction
    between 5-fluorocyclopenta-1,3-diene and maleic anhydride.
    """
    question = "5-fluorocyclopenta-1,3-diene is reacted with maleic anhydride. What is the major product?"
    options = {
        "A": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione"
    }
    llm_answer = "C"

    # --- Step 1: Endo vs. Exo Selectivity ---
    # The Alder Endo Rule states that for kinetically controlled Diels-Alder reactions,
    # the endo product is favored due to stabilizing secondary orbital interactions.
    # We need to map the IUPAC names to endo/exo structures.
    # For the bicyclo[2.2.1]heptane system fused with the anhydride:
    # Endo skeleton: (3aR,4S,7R,7aS) or its enantiomer.
    # Exo skeleton: (3aR,4R,7S,7aS) or its enantiomer.
    
    endo_pattern = r"\(3aR,4S,7R,7aS,.*\)"
    exo_pattern = r"\(3aR,4R,7S,7aS,.*\)"

    possible_options = {}
    for option, name in options.items():
        if re.search(endo_pattern, name):
            possible_options[option] = "endo"
        elif re.search(exo_pattern, name):
            possible_options[option] = "exo"

    # The major product should be endo.
    major_product_orientation = "endo"
    
    # Filter out exo products
    endo_options = {k: v for k, v in possible_options.items() if v == major_product_orientation}
    
    if not endo_options:
        return "Error in analysis: No options correspond to the expected endo product."
    if llm_answer not in endo_options:
        return f"Incorrect. The answer {llm_answer} corresponds to an exo product. The Diels-Alder reaction under kinetic control favors the endo product due to the Alder Endo Rule."

    # --- Step 2: Facial Selectivity (Syn vs. Anti Attack) ---
    # For 5-substituted cyclopentadienes, small, electronegative substituents like Fluorine
    # favor syn-facial attack due to electronic effects (e.g., negative hyperconjugation)
    # that outweigh steric hindrance.
    favored_attack = "syn-facial"

    # --- Step 3: Translate Attack Pathway to Product Structure ---
    # It's crucial to map the attack pathway to the final product's relative stereochemistry.
    # - syn-facial attack -> anti-adduct (Fluorine is anti to the anhydride ring)
    # - anti-facial attack -> syn-adduct (Fluorine is syn to the anhydride ring)
    major_adduct_type = "anti-adduct"

    # --- Step 4: Map Product Structure to IUPAC Nomenclature (8r/8s) ---
    # We need to determine if the 'anti-adduct' corresponds to '8r' or '8s'.
    # This requires a Cahn-Ingold-Prelog (CIP) analysis for the C8 pseudoasymmetric center.
    # Priorities on C8: (1)F > (2)C7(R) > (3)C4(S) > (4)H.
    # In the anti-adduct, F is on the opposite side of the anhydride ring. Viewing down the
    # C8-H bond (lowest priority away), the sequence 1->2->3 is clockwise.
    # Clockwise = R configuration.
    # Therefore, the anti-adduct corresponds to the '8r' descriptor.
    
    anti_adduct_descriptor = "8r"
    
    # Find the option that is both 'endo' and has the '8r' descriptor.
    correct_option = None
    for option_key in endo_options:
        if anti_adduct_descriptor in options[option_key]:
            correct_option = option_key
            break
            
    if correct_option is None:
        return "Error in analysis: No endo option matches the derived '8r' descriptor for the anti-adduct."

    # --- Step 5: Final Check ---
    if llm_answer == correct_option:
        return "Correct"
    else:
        reason = f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option}.\n"
        reason += f"Reasoning:\n"
        reason += f"1. Endo/Exo Selectivity: The reaction favors the 'endo' product, eliminating options A and D.\n"
        reason += f"2. Facial Selectivity: Electronic effects favor 'syn-facial' attack for the fluorine substituent.\n"
        reason += f"3. Product Structure: Syn-facial attack leads to the 'anti-adduct' (fluorine is anti to the anhydride ring).\n"
        reason += f"4. Nomenclature: CIP analysis shows the 'anti-adduct' corresponds to the '8r' descriptor.\n"
        reason += f"Therefore, the major product is the endo, anti-adduct, which is option {correct_option}."
        return reason

# Run the check
result = check_diels_alder_product()
print(result)