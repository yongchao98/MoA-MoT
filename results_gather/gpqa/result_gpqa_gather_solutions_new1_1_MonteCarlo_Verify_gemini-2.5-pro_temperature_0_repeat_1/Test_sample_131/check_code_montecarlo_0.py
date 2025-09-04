import collections

def check_nmr_answer():
    """
    Checks the correctness of the answer to the NMR spectroscopy question.

    The function defines the predicted 1H NMR spectra for four candidate compounds.
    It then simulates the spectrum for a 1:1 mixture for each option and compares
    it against the experimental data given in the question.
    """

    # --- Step 1: Define predicted NMR data for each compound ---
    # Based on chemical structure and symmetry.
    # 's' = singlet, 'd' = doublet, 't' = triplet, 'q' = quartet
    # Integrals are relative proton counts.
    # Note: The analysis for 1,2,3,5-tetramethylbenzene (Isodurene) correctly identifies
    # two non-equivalent aromatic protons, leading to two 1H singlets.
    # Note: The two adjacent aromatic protons in 1,2,3,4-tetramethylbenzene (Prehnitene)
    # form an AA' system that is often observed as a singlet, which is consistent
    # with the problem's simplification.
    compounds_data = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"multiplicity": ['s'], "integrals": [2]},
            "aliphatic": {"multiplicity": ['s'], "integrals": [12]}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"multiplicity": ['s', 's'], "integrals": [1, 1]},
            "aliphatic": {"multiplicity": ['s', 's', 's'], "integrals": [6, 3, 3]}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"multiplicity": ['s'], "integrals": [2]},
            "aliphatic": {"multiplicity": ['s', 's'], "integrals": [6, 6]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"multiplicity": ['s'], "integrals": [4]},
            "aliphatic": {"multiplicity": ['q', 't'], "integrals": [4, 6]}
        }
    }

    # --- Step 2: Define the options and the target spectrum from the question ---
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"] # Note: The original prompt had a different option D, but this is the one from the LLM answers. We will check all combinations.
    }
    
    # The question prompt's option D is 1,2,3,5-tetramethylbenzene and 1,4-diethylbenzene.
    # Let's ensure the options dictionary matches the question prompt.
    options["D"] = ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"]


    target_spectrum = {
        "aromatic": {"signals": 2, "multiplicity": "singlets", "ratio": [1, 1]},
        "aliphatic": {"signals": 3, "multiplicity": "singlets", "ratio": [2, 1, 1]}
    }
    
    llm_answer = "A"

    # --- Step 3: Function to check if a ratio matches ---
    def check_ratio(integrals, target_ratio):
        if len(integrals) != len(target_ratio):
            return False
        
        sorted_integrals = sorted(integrals)
        sorted_target = sorted(target_ratio)
        
        # Ratios must be proportional
        factor = sorted_integrals[0] / sorted_target[0]
        if factor == 0: return False # Avoid division by zero if integrals are 0
        
        for i in range(len(sorted_integrals)):
            if abs(sorted_integrals[i] / sorted_target[i] - factor) > 1e-9:
                return False
        return True

    # --- Step 4: Evaluate each option ---
    correct_option = None
    for option_key, compound_names in options.items():
        c1_name, c2_name = compound_names
        c1_data = compounds_data[c1_name]
        c2_data = compounds_data[c2_name]

        # Combine spectra for the 1:1 mixture
        
        # Check Aliphatic Multiplicity
        combined_aliphatic_mult = c1_data["aliphatic"]["multiplicity"] + c2_data["aliphatic"]["multiplicity"]
        if not all(m == 's' for m in combined_aliphatic_mult):
            continue # This option is invalid as it has non-singlet aliphatic signals

        # Check Aromatic Signals
        combined_aromatic_integrals = c1_data["aromatic"]["integrals"] + c2_data["aromatic"]["integrals"]
        if len(combined_aromatic_integrals) != target_spectrum["aromatic"]["signals"]:
            continue
        if not check_ratio(combined_aromatic_integrals, target_spectrum["aromatic"]["ratio"]):
            continue

        # Check Aliphatic Signals
        combined_aliphatic_integrals = c1_data["aliphatic"]["integrals"] + c2_data["aliphatic"]["integrals"]
        if len(combined_aliphatic_integrals) != target_spectrum["aliphatic"]["signals"]:
            continue
        if not check_ratio(combined_aliphatic_integrals, target_spectrum["aliphatic"]["ratio"]):
            continue
            
        # If all checks pass, this is the correct option
        correct_option = option_key
        break # Found the match

    # --- Step 5: Final Verification ---
    if correct_option == llm_answer:
        return "Correct"
    elif correct_option is None:
        return f"The provided answer '{llm_answer}' is incorrect. My analysis could not find any option that perfectly matches the criteria without assuming signal overlap. The closest match is Option A, but there is a discrepancy in the analysis of some LLMs. The provided answer A is the only one that works without assuming accidental degeneracy."
    else:
        return f"The provided answer '{llm_answer}' is incorrect. The correct answer is '{correct_option}'. Reason: Option {llm_answer} does not produce the correct number of signals or integration ratios. Option {correct_option} is the only one that perfectly matches the experimental data: two aromatic singlets in a 1:1 ratio (from two compounds each giving a 2H singlet) and three aliphatic singlets in a 2:1:1 ratio (from one compound giving a 12H singlet and the other giving two 6H singlets)."

# Execute the check
result = check_nmr_answer()
print(result)