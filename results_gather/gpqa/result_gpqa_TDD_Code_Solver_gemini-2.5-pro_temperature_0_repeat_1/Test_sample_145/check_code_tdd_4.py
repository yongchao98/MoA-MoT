import re

def check_diels_alder_stereochemistry():
    """
    This function checks the correctness of the provided answer for the Diels-Alder
    reaction between 5-fluorocyclopenta-1,3-diene and maleic anhydride.

    The verification is based on two key stereochemical principles:
    1.  **Endo/Exo Selectivity**: The Diels-Alder reaction is under kinetic control
        and strongly favors the 'endo' adduct due to stabilizing secondary orbital
        overlap between the diene and dienophile.
    2.  **Facial (Syn/Anti) Selectivity**: The dienophile (maleic anhydride) can
        attack the diene from the same side ('syn') as the fluorine substituent or
        from the opposite side ('anti'). While fluorine can have a syn-directing
        electronic effect, the 'endo-syn' transition state would suffer from
        severe steric repulsion between the fluorine atom and the anhydride ring.
        To avoid this prohibitive steric clash, the reaction proceeds via 'anti'
        attack.

    Therefore, the major product is the 'endo-anti' adduct.

    The function maps these chemical terms to the provided IUPAC names:
    -   **Endo Adduct**: Contains the stereodescriptors (...,4R,7S,...).
    -   **Exo Adduct**: Contains the stereodescriptors (...,4S,7R,...).
    -   **Anti Product** (Fluorine is anti to the anhydride bridge): Contains (...,8r).
    -   **Syn Product** (Fluorine is syn to the anhydride bridge): Contains (...,8s).
    """
    llm_answer_choice = "B"

    options = {
        "A": "(3aR,4S,7R,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "B": "(3aR,4R,7S,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "C": "(3aR,4R,7S,7aS,8s)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
        "D": "(3aR,4S,7R,7aS,8r)-8-fluoro-3a,4,7,7a-tetrahydro-4,7-methanoisobenzofuran-1,3-dione",
    }

    # Step 1: Identify candidates based on the 'endo' rule.
    # The endo product has the (4R,7S) configuration.
    endo_candidates = {key for key, name in options.items() if "4R,7S" in name}

    # Step 2: Identify candidates based on the 'anti' attack rule.
    # The anti product has the (8r) configuration.
    anti_candidates = {key for key, name in options.items() if re.search(r',8r\)', name)}

    # Step 3: The major product is the intersection of both favored conditions (endo and anti).
    major_product_keys = endo_candidates.intersection(anti_candidates)

    # There should be only one product that satisfies both conditions.
    if len(major_product_keys) != 1:
        return (f"Constraint check failed: The chemical principles did not lead to a unique answer. "
                f"Endo candidates: {endo_candidates}. Anti candidates: {anti_candidates}. "
                f"The intersection is {major_product_keys}, which does not contain a single option.")

    correct_choice = major_product_keys.pop()

    # Step 4: Compare the determined correct choice with the LLM's answer.
    if llm_answer_choice == correct_choice:
        return "Correct"
    else:
        # Determine the stereochemistry of the incorrect answer provided by the LLM.
        llm_choice_name = options[llm_answer_choice]
        is_endo = "4R,7S" in llm_choice_name
        is_anti = re.search(r',8r\)', llm_choice_name) is not None
        
        llm_product_type = f"{'endo' if is_endo else 'exo'}-{'anti' if is_anti else 'syn'}"

        return (f"Incorrect. The major product of this reaction is the 'endo-anti' adduct. "
                f"This corresponds to option {correct_choice}. The provided answer, {llm_answer_choice}, "
                f"is the '{llm_product_type}' adduct. This is incorrect because the reaction favors the 'endo' "
                f"pathway for kinetic reasons and the 'anti' attack to avoid severe steric hindrance.")

# Run the check
result = check_diels_alder_stereochemistry()
print(result)