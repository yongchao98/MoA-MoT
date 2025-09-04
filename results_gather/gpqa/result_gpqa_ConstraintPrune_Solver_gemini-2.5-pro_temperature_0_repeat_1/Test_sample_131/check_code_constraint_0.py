import collections

def check_nmr_answer():
    """
    This function verifies the correct pair of compounds based on 1H NMR data.
    It models the properties of each potential compound and checks which pair
    matches the spectral constraints given in the problem.
    """

    # Step 1: Define the predicted 1H NMR data for each candidate compound.
    # Data format: {name: {'formula': str, 'aromatic': [(multiplicity, integration), ...], 'aliphatic': [...]}}
    # This data is derived from the chemical structure and symmetry of each molecule.
    candidates = {
        "1,2,4,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 2)],  # Two equivalent aromatic H's -> 1 singlet (2H)
            "aliphatic": [("singlet", 12)] # Four equivalent methyl groups -> 1 singlet (12H)
        },
        "1,2,3,5-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 1), ("singlet", 1)], # Two non-equivalent aromatic H's -> 2 singlets (1H each)
            "aliphatic": [("singlet", 6), ("singlet", 3), ("singlet", 3)] # Three non-equivalent methyl groups -> 3 singlets (6H, 3H, 3H)
        },
        "1,2,3,4-tetramethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 2)], # Two equivalent aromatic H's -> 1 singlet (2H)
            "aliphatic": [("singlet", 6), ("singlet", 6)] # Two sets of equivalent methyl groups -> 2 singlets (6H each)
        },
        "1,4-diethylbenzene": {
            "formula": "C10H14",
            "aromatic": [("singlet", 4)], # Four equivalent aromatic H's -> 1 singlet (4H)
            "aliphatic": [("quartet", 4), ("triplet", 6)] # CH2 and CH3 of ethyl groups -> 1 quartet (4H), 1 triplet (6H)
        }
    }

    # Step 2: Define the options and the proposed answer from the LLM.
    options = {
        "A": ("1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"),
        "B": ("1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"),
        "C": ("1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"),
        "D": ("1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene")
    }
    
    llm_answer_key = "B"
    
    # Step 3: Retrieve the compounds for the given answer.
    try:
        compound1_name, compound2_name = options[llm_answer_key]
        c1_data = candidates[compound1_name]
        c2_data = candidates[compound2_name]
    except KeyError:
        return f"Invalid answer key '{llm_answer_key}'. It is not one of the options A, B, C, or D."

    # Step 4: Check the constraints for the selected pair.
    
    # Constraint 1: Molecular Formula must be C10H14 for both.
    if not (c1_data["formula"] == "C10H14" and c2_data["formula"] == "C10H14"):
        return "Incorrect. Constraint violated: Both compounds must have the molecular formula C10H14."

    # Combine the signals for a 1:1 mixture.
    combined_aromatic = c1_data["aromatic"] + c2_data["aromatic"]
    combined_aliphatic = c1_data["aliphatic"] + c2_data["aliphatic"]

    # Constraint 2: Aromatic region must have two singlets in a 1:1 ratio.
    if len(combined_aromatic) != 2:
        return f"Incorrect. Aromatic region constraint violated: The mixture should have 2 signals, but the proposed pair has {len(combined_aromatic)}."
    if not all(sig[0] == "singlet" for sig in combined_aromatic):
        return "Incorrect. Aromatic region constraint violated: All signals in this region must be singlets."
    # For a 1:1 ratio, the two integrations must be equal.
    if combined_aromatic[0][1] != combined_aromatic[1][1]:
        return f"Incorrect. Aromatic region constraint violated: The signals should have a 1:1 integration ratio, but the proposed pair gives a ratio of {combined_aromatic[0][1]}:{combined_aromatic[1][1]}."

    # Constraint 3: Aliphatic region must have three singlets in a 2:1:1 ratio.
    if not all(sig[0] == "singlet" for sig in combined_aliphatic):
        # This check is important for eliminating 1,4-diethylbenzene.
        non_singlets = [sig[0] for sig in combined_aliphatic if sig[0] != 'singlet']
        return f"Incorrect. Aliphatic region constraint violated: All signals must be singlets, but the proposed pair contains a {non_singlets[0]}."
    if len(combined_aliphatic) != 3:
        return f"Incorrect. Aliphatic region constraint violated: The mixture should have 3 signals, but the proposed pair has {len(combined_aliphatic)}."
    
    # Check the 2:1:1 ratio.
    aliphatic_integrations = sorted([sig[1] for sig in combined_aliphatic])
    unit = aliphatic_integrations[0] # The smallest integration should be the '1' part of the ratio.
    expected_ratio = sorted([2 * unit, 1 * unit, 1 * unit])
    if aliphatic_integrations != expected_ratio:
        return f"Incorrect. Aliphatic region constraint violated: The signals should have a 2:1:1 integration ratio, but the proposed pair gives integrations of {aliphatic_integrations}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the checker and print the result.
result = check_nmr_answer()
print(result)