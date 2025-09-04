def check_correctness():
    """
    This function checks the correctness of the final answer by systematically applying
    the spectral data constraints to each candidate molecule.
    """

    # Define the properties of each candidate molecule based on its structure.
    candidates = {
        "A": {
            "name": "Phenyl chloroformate",
            "has_carboxylic_acid": False,
            "aromatic_protons": 5,
            "substitution_pattern": "monosubstituted"
        },
        "B": {
            "name": "2-chlorobenzoic acid",
            "has_carboxylic_acid": True,
            "aromatic_protons": 4,
            "substitution_pattern": "ortho"  # 1,2-disubstituted
        },
        "C": {
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "has_carboxylic_acid": False, # It has phenol and aldehyde groups
            "aromatic_protons": 3,
            "substitution_pattern": "trisubstituted"
        },
        "D": {
            "name": "4-chlorobenzoic acid",
            "has_carboxylic_acid": True,
            "aromatic_protons": 4,
            "substitution_pattern": "para"  # 1,4-disubstituted
        }
    }

    # The final answer provided by the LLM to be checked.
    llm_answer_key = "D"

    # --- Define Constraints from Spectral Data ---
    # MS: MW ~156 with 1 Cl atom. All candidates have formula C7H5ClO2, so all pass this.
    # IR: Broad peak 3500-2700 cm-1 + sharp peak 1720 cm-1 => Carboxylic Acid
    ir_requires_carboxylic_acid = True
    # NMR: 11.0 ppm (s, 1H) => Carboxylic Acid
    # NMR: 8.02 ppm (d, 2H) + 7.72 (d, 2H) => 4 aromatic protons, para-substitution pattern
    nmr_requires_para_substitution = True
    nmr_requires_4_aromatic_protons = True

    # --- Verification Logic ---
    passing_candidates = []
    failure_reasons = {}

    for key, molecule in candidates.items():
        # Check IR constraint
        ir_ok = ir_requires_carboxylic_acid == molecule["has_carboxylic_acid"]
        
        # Check NMR constraints
        nmr_pattern_ok = (nmr_requires_para_substitution and molecule["substitution_pattern"] == "para")
        nmr_proton_count_ok = (nmr_requires_4_aromatic_protons and molecule["aromatic_protons"] == 4)
        # The carboxylic acid proton signal is redundant with the IR check
        nmr_ok = nmr_pattern_ok and nmr_proton_count_ok

        if ir_ok and nmr_ok:
            passing_candidates.append(key)
        else:
            reasons = []
            if not ir_ok:
                reasons.append("fails IR constraint (does not have a carboxylic acid group)")
            if not nmr_proton_count_ok:
                reasons.append(f"fails NMR proton count (has {molecule['aromatic_protons']} aromatic protons, expected 4)")
            if not nmr_pattern_ok:
                reasons.append(f"fails NMR substitution pattern (is {molecule['substitution_pattern']}, expected para)")
            failure_reasons[key] = "; ".join(reasons)

    # --- Final Judgement ---
    if llm_answer_key in passing_candidates and len(passing_candidates) == 1:
        return "Correct"
    elif llm_answer_key not in passing_candidates:
        reason = f"Incorrect. The proposed answer '{llm_answer_key}' ({candidates[llm_answer_key]['name']}) is wrong because it {failure_reasons.get(llm_answer_key, 'does not fit the data')}.\n"
        if passing_candidates:
            correct_key = passing_candidates[0]
            reason += f"The only molecule that satisfies all constraints is '{correct_key}' ({candidates[correct_key]['name']})."
        else:
            reason += "No single candidate satisfies all the given constraints."
        return reason
    else: # This case handles ambiguity where the LLM answer is one of multiple correct options
        return f"Incorrect. The answer '{llm_answer_key}' is plausible, but the data is ambiguous as candidates {passing_candidates} all satisfy the constraints."

# The final output of the check.
result = check_correctness()
print(result)