def check_spectroscopy_answer():
    """
    This function checks the correctness of the identified chemical structure based on spectral data.
    It defines the properties of four candidate molecules and verifies them against the given
    Mass Spec, IR, and NMR data.
    """

    # --- Define Spectral Data Constraints ---
    # MS: M+ at 156, M+2 at 158 (~3:1 ratio) -> MW=156, contains one Cl
    # IR: Broad 3500-2700 cm-1 & sharp 1720 cm-1 -> Carboxylic Acid
    # NMR: 11.0 ppm (1H) -> Carboxylic Acid proton
    # NMR: Two doublets, each 2H -> 4 aromatic protons, para-disubstituted pattern

    # --- Define Candidate Molecules ---
    candidates = {
        "A": {
            "name": "2-chlorobenzoic acid",
            "mw_35cl": 156,
            "has_cl": True,
            "is_carboxylic_acid": True,
            "substitution_pattern": "ortho", # 1,2-disubstituted
            "aromatic_protons": 4
        },
        "B": {
            "name": "Phenyl chloroformate",
            "mw_35cl": 156,
            "has_cl": True,
            "is_carboxylic_acid": False, # It's a chloroformate
            "substitution_pattern": "monosubstituted",
            "aromatic_protons": 5
        },
        "C": {
            "name": "4-chlorobenzoic acid",
            "mw_35cl": 156,
            "has_cl": True,
            "is_carboxylic_acid": True,
            "substitution_pattern": "para", # 1,4-disubstituted
            "aromatic_protons": 4
        },
        "D": {
            "name": "3-Chloro-2-hydroxybenzaldehyde",
            "mw_35cl": 156,
            "has_cl": True,
            "is_carboxylic_acid": False, # It's an aldehyde and a phenol
            "substitution_pattern": "trisubstituted",
            "aromatic_protons": 3
        }
    }

    # The final answer from the LLM to be checked
    llm_answer_choice = "C"
    
    chosen_candidate = candidates.get(llm_answer_choice)

    # --- Verification Steps ---

    # 1. Mass Spec Check: All candidates have MW 156 and one Cl, so this is a basic check.
    if chosen_candidate["mw_35cl"] != 156 or not chosen_candidate["has_cl"]:
        return f"Incorrect. The MS data requires a molecule with MW=156 containing one chlorine atom. {chosen_candidate['name']} does not match this."

    # 2. Functional Group Check (IR & NMR acid proton)
    # The IR data (broad 3500-2700 cm-1) and NMR data (11.0 ppm signal) are definitive for a carboxylic acid.
    if not chosen_candidate["is_carboxylic_acid"]:
        return f"Incorrect. The IR and NMR data strongly indicate a carboxylic acid functional group. {chosen_candidate['name']} is not a carboxylic acid."

    # 3. Aromatic Proton Count Check (NMR integration)
    # The NMR data shows a total of 4 aromatic protons (2H + 2H).
    if chosen_candidate["aromatic_protons"] != 4:
        return f"Incorrect. The NMR data shows 4 aromatic protons. {chosen_candidate['name']} has {chosen_candidate['aromatic_protons']} aromatic protons."

    # 4. Substitution Pattern Check (NMR splitting)
    # The NMR aromatic signals (two doublets, each for 2H) indicate a symmetrical para (1,4) substitution.
    if chosen_candidate["substitution_pattern"] != "para":
        return f"Incorrect. The NMR splitting pattern (two doublets for 2H each) is characteristic of a para-substituted ring. {chosen_candidate['name']} has a(n) '{chosen_candidate['substitution_pattern']}' substitution pattern, which would produce a different spectrum."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_spectroscopy_answer()
print(result)