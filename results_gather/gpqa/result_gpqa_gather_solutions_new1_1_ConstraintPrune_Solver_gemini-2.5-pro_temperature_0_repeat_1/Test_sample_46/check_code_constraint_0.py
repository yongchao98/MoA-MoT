def check_spectroscopy_answer():
    """
    This function checks the correctness of the identified compound based on spectral data.
    It simulates the process of structure elucidation by applying constraints from IR and NMR data
    to a list of candidate molecules.
    """

    # --- 1. Define Spectral Data Constraints ---
    # IR Data Interpretation
    # Two N-H bands (3420, 3325 cm-1) strongly indicate a primary amine (-NH2).
    ir_constraint_primary_amine = True
    # C=O stretch at 1720 cm-1 is characteristic of a conjugated ester. Amides are typically lower (~1650-1690 cm-1).
    ir_constraint_ester = True

    # NMR Data Interpretation
    # Two doublets in the aromatic region (7.0, 8.0 ppm), each 2H, indicate 1,4- (para) substitution.
    nmr_constraint_para_substitution = True
    # A triplet (~1.2 ppm) and a quartet (~4.5 ppm) indicate an ethyl group.
    # The high chemical shift of the quartet (4.5 ppm) is specific to an ethyl group attached to an ester oxygen (-O-CH2CH3).
    nmr_constraint_ethyl_ester = True
    # A broad singlet at 4.0 ppm (2H) confirms the primary amine.
    nmr_constraint_primary_amine = True

    # --- 2. Define Candidate Molecules and Their Properties ---
    candidates = {
        "A": {
            "name": "ethyl 4-aminobenzoate",
            "has_primary_amine": True,
            "is_ester": True,
            "is_amide": False,
            "substitution": "para",
            "has_ethyl_ester_group": True, # The key feature: -COOCH2CH3
        },
        "B": {
            "name": "4-aminophenyl propionate",
            "has_primary_amine": True,
            "is_ester": True,
            "is_amide": False,
            "substitution": "para",
            # Has a propionate group (-OCOCH2CH3). The -CH2- quartet would be ~2.5 ppm, not 4.5 ppm.
            "has_ethyl_ester_group": False,
        },
        "C": {
            "name": "3-ethoxybenzamide",
            "has_primary_amine": False, # It's a primary amide, not a primary amine on the ring.
            "is_ester": False,
            "is_amide": True,
            "substitution": "meta", # 1,3-substitution
            "has_ethyl_ester_group": False,
        },
        "D": {
            "name": "N-(4-ethoxyphenyl)formamide",
            "has_primary_amine": False, # It's a secondary amide.
            "is_ester": False,
            "is_amide": True,
            "substitution": "para",
            "has_ethyl_ester_group": False,
        }
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = "A"

    # --- 3. Evaluate Each Candidate Against Constraints ---
    valid_candidates = []
    reasons_for_elimination = {}

    for key, properties in candidates.items():
        errors = []
        # Check IR constraints
        if ir_constraint_primary_amine and not properties["has_primary_amine"]:
            errors.append("Fails IR constraint: Does not have a primary amine group (-NH2).")
        if ir_constraint_ester and not properties["is_ester"]:
            errors.append("Fails IR constraint: Is not an ester (C=O at 1720 cm-1).")
        
        # Check NMR constraints
        if nmr_constraint_para_substitution and properties["substitution"] != "para":
            errors.append("Fails NMR constraint: Is not para-substituted (aromatic region is not two doublets).")
        if nmr_constraint_ethyl_ester and not properties["has_ethyl_ester_group"]:
            errors.append("Fails NMR constraint: Does not have an ethyl ester group (-O-CH2CH3) which would cause a quartet at ~4.5 ppm.")

        if not errors:
            valid_candidates.append(key)
        else:
            reasons_for_elimination[key] = " ".join(errors)

    # --- 4. Determine Correctness of the LLM's Answer ---
    if len(valid_candidates) == 1:
        correct_answer = valid_candidates[0]
        if llm_answer == correct_answer:
            return "Correct"
        else:
            return f"Incorrect. The provided answer was '{llm_answer}', but the only structure that fits all spectral data is '{correct_answer}' ({candidates[correct_answer]['name']})."
    elif len(valid_candidates) == 0:
        return f"Incorrect. No candidate fits all the data. The provided answer was '{llm_answer}'. Reasons for elimination: {reasons_for_elimination}"
    else:
        return f"Incorrect. The data is ambiguous as multiple candidates {valid_candidates} fit the criteria. The provided answer was '{llm_answer}'."

# Run the check and print the result
result = check_spectroscopy_answer()
print(result)