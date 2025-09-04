def check_answer_correctness():
    """
    This function checks the correctness of the given answer by verifying if the proposed molecule's
    properties match the constraints derived from the spectral data.
    """
    # The final answer provided by the LLM
    llm_answer = "B"

    # Database of molecular properties for each option
    molecules = {
        "A": {
            "name": "Phenyl chloroformate",
            "mw": 156.57,
            "has_cl": True,
            "has_carboxylic_acid": False,
            "aromatic_pattern": "monosubstituted"
        },
        "B": {
            "name": "4-chlorobenzoic acid",
            "mw": 156.57,
            "has_cl": True,
            "has_carboxylic_acid": True,
            "aromatic_pattern": "para"
        },
        "C": {
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "mw": 156.57,
            "has_cl": True,
            "has_carboxylic_acid": False, # Has phenol and aldehyde, not carboxylic acid
            "aromatic_pattern": "trisubstituted"
        },
        "D": {
            "name": "2-chlorobenzoic acid",
            "mw": 156.57,
            "has_cl": True,
            "has_carboxylic_acid": True,
            "aromatic_pattern": "ortho"
        }
    }

    # Constraints derived from the spectral data in the question
    constraints = {
        "mw_peak": 156,
        "has_cl": True,
        "has_carboxylic_acid": True, # From IR (broad 3500-2700, 1720) and NMR (11.0 ppm)
        "aromatic_pattern": "para" # From NMR (two doublets, 2H each)
    }

    # Select the candidate molecule based on the LLM's answer
    candidate = molecules.get(llm_answer)

    if not candidate:
        return f"Invalid answer key '{llm_answer}'. Please choose from {list(molecules.keys())}."

    # --- Perform Checks ---
    errors = []

    # Check 1: Mass Spec (MW and Chlorine presence)
    # All options have the correct MW and a chlorine atom, so this check is implicitly passed.
    if not (constraints["mw_peak"] - 1 < candidate["mw"] < constraints["mw_peak"] + 1):
        errors.append(f"MS constraint failed: Expected MW around {constraints['mw_peak']}, but {candidate['name']} has MW {candidate['mw']}.")
    if constraints["has_cl"] and not candidate["has_cl"]:
        errors.append("MS constraint failed: The M+2 peak indicates a Chlorine atom, which is missing.")

    # Check 2: IR and NMR Functional Group (Carboxylic Acid)
    if constraints["has_carboxylic_acid"] and not candidate["has_carboxylic_acid"]:
        errors.append(f"Functional group constraint failed: The IR and NMR data strongly indicate a carboxylic acid, which is absent in {candidate['name']}.")

    # Check 3: NMR Aromatic Pattern
    if candidate["aromatic_pattern"] != constraints["aromatic_pattern"]:
        errors.append(f"NMR aromatic pattern constraint failed: The NMR shows a pattern for a '{constraints['aromatic_pattern']}' substituted ring, but {candidate['name']} has a '{candidate['aromatic_pattern']}' substitution pattern.")

    # --- Final Verdict ---
    if not errors:
        return "Correct"
    else:
        # Join all found errors into a single reason string
        reason = "Incorrect. The proposed answer does not satisfy the following constraints:\n"
        for i, error in enumerate(errors):
            reason += f"- {error}\n"
        return reason.strip()

# Run the check and print the result
result = check_answer_correctness()
print(result)