def check_answer_correctness():
    """
    This function checks the correctness of the final answer by verifying
    the properties of the proposed molecule against the given spectral data.
    """
    # --- Spectral Data Constraints from the Question ---
    constraints = {
        "molecular_weight": 156,
        "has_one_chlorine": True,
        "is_carboxylic_acid": True,  # From IR (broad 3500-2700, C=O 1720) and NMR (11.0 ppm)
        "aromatic_protons_count": 4, # From NMR (2H + 2H)
        "substitution_pattern": "para" # From NMR (two doublets for 2H each)
    }

    # --- Properties of the Candidate Molecules ---
    molecules = {
        "A": {
            "name": "Phenyl chloroformate",
            "mw": 156,
            "has_cl": True,
            "is_carboxylic_acid": False,
            "aromatic_protons": 5,
            "substitution": "monosubstituted"
        },
        "B": {
            "name": "4-chlorobenzoic acid",
            "mw": 156,
            "has_cl": True,
            "is_carboxylic_acid": True,
            "aromatic_protons": 4,
            "substitution": "para"
        },
        "C": {
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "mw": 156,
            "has_cl": True,
            "is_carboxylic_acid": False,
            "aromatic_protons": 3,
            "substitution": "trisubstituted"
        },
        "D": {
            "name": "2-chlorobenzoic acid",
            "mw": 156,
            "has_cl": True,
            "is_carboxylic_acid": True,
            "aromatic_protons": 4,
            "substitution": "ortho"
        }
    }

    # The final answer provided by the LLM analysis is 'B'.
    proposed_answer_key = "B"
    molecule_to_check = molecules[proposed_answer_key]

    # --- Verification Steps ---

    # 1. Check for Carboxylic Acid functional group
    if molecule_to_check["is_carboxylic_acid"] != constraints["is_carboxylic_acid"]:
        return (f"Incorrect. The IR and NMR data confirm the presence of a carboxylic acid, "
                f"but {molecule_to_check['name']} is not a carboxylic acid.")

    # 2. Check the number of aromatic protons
    if molecule_to_check["aromatic_protons"] != constraints["aromatic_protons_count"]:
        return (f"Incorrect. The 1H NMR spectrum shows {constraints['aromatic_protons_count']} aromatic protons (2H + 2H), "
                f"but {molecule_to_check['name']} has {molecule_to_check['aromatic_protons']} aromatic protons.")

    # 3. Check the substitution pattern (the most discriminating factor)
    if molecule_to_check["substitution"] != constraints["substitution_pattern"]:
        return (f"Incorrect. The 1H NMR spectrum shows a symmetrical pattern of two doublets, "
                f"which is characteristic of a '{constraints['substitution_pattern']}' substituted ring. "
                f"However, {molecule_to_check['name']} has an '{molecule_to_check['substitution']}' substitution pattern, "
                f"which would produce a different NMR spectrum.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_answer_correctness()
print(result)