def check_answer_correctness():
    """
    This function checks the correctness of the final answer by systematically applying
    the spectral data constraints to the candidate molecules.
    """

    # Define the properties of the candidate molecules based on chemical knowledge.
    candidates = [
        {
            "id": "A",
            "name": "2-chlorobenzoic acid",
            "has_chlorine": True,
            "nominal_mass": 156,
            "functional_group": "Carboxylic Acid",
            "aromatic_pattern": "ortho"  # 1,2-disubstituted, would give 4 complex signals for 4H
        },
        {
            "id": "B",
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "has_chlorine": True,
            "nominal_mass": 156,
            "functional_group": "Aldehyde/Phenol", # Not a carboxylic acid
            "aromatic_pattern": "trisubstituted" # Would give 3 signals for 3H
        },
        {
            "id": "C",
            "name": "Phenyl chloroformate",
            "has_chlorine": True,
            "nominal_mass": 156,
            "functional_group": "Chloroformate", # Not a carboxylic acid, no O-H group
            "aromatic_pattern": "monosubstituted" # Would give signals for 5H
        },
        {
            "id": "D",
            "name": "4-chlorobenzoic acid",
            "has_chlorine": True,
            "nominal_mass": 156,
            "functional_group": "Carboxylic Acid",
            "aromatic_pattern": "para" # 1,4-disubstituted, gives 2 doublets for 4H
        }
    ]

    # The final answer provided by the LLM to be checked.
    llm_final_answer_id = "D"

    # --- Constraint 1: Mass Spectrometry ---
    # Data: M+ at m/z=156, M+2 at m/z=158 (32% intensity).
    # Interpretation: The molecule has a nominal mass of 156 and contains one chlorine atom.
    survivors = [c for c in candidates if c["has_chlorine"] and c["nominal_mass"] == 156]
    # In this case, all candidates pass this check.

    # --- Constraint 2: IR Spectroscopy ---
    # Data: Broad peak 3500-2700 cm⁻¹, strong sharp peak at 1720 cm⁻¹.
    # Interpretation: The molecule contains a carboxylic acid group.
    ir_survivors = []
    for c in survivors:
        if c["functional_group"] == "Carboxylic Acid":
            ir_survivors.append(c)
        else:
            # Check if the LLM's answer was eliminated here
            if c["id"] == llm_final_answer_id:
                return f"Incorrect. The provided answer {c['id']} ({c['name']}) is wrong because it is not a carboxylic acid, which is strongly indicated by the IR spectrum (broad O-H stretch and C=O peak at 1720 cm⁻¹)."
    survivors = ir_survivors

    # --- Constraint 3: 1H NMR Spectroscopy ---
    # Data: 11.0 ppm (s, 1H) -> Confirms carboxylic acid.
    # Data: 8.02 ppm (d, 2H), 7.72 (d, 2H) -> 4 aromatic protons in a para-substituted pattern.
    # Interpretation: The aromatic ring is 1,4-disubstituted (para).
    nmr_survivors = []
    for c in survivors:
        if c["aromatic_pattern"] == "para":
            nmr_survivors.append(c)
        else:
            # Check if the LLM's answer was eliminated here
            if c["id"] == llm_final_answer_id:
                return f"Incorrect. The provided answer {c['id']} ({c['name']}) is wrong because its '{c['aromatic_pattern']}' substitution pattern does not match the ¹H NMR data, which shows a classic 'para' pattern (two doublets, each for 2H)."
    survivors = nmr_survivors

    # --- Final Verification ---
    if len(survivors) == 1:
        correct_candidate = survivors[0]
        if correct_candidate["id"] == llm_final_answer_id:
            return "Correct"
        else:
            return f"Incorrect. The analysis points to {correct_candidate['id']} ({correct_candidate['name']}), but the provided answer was {llm_final_answer_id}."
    elif len(survivors) == 0:
        return "Incorrect. No candidate satisfies all the given spectral constraints."
    else:
        surviving_names = [c['name'] for c in survivors]
        return f"Incorrect. The data is ambiguous as multiple candidates remain: {', '.join(surviving_names)}."

# Run the check
result = check_answer_correctness()
print(result)