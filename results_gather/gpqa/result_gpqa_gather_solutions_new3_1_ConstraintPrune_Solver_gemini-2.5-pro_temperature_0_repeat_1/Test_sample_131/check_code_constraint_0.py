import collections

def check_nmr_answer():
    """
    Checks the correctness of the answer to the NMR problem by simulating the spectra of the mixtures.
    """
    # Define the predicted 1H NMR spectral properties for each compound based on its structure and symmetry.
    # Note: The analysis assumes that the two adjacent aromatic protons in 1,2,3,4-tetramethylbenzene
    # produce a singlet, a common simplification in such problems.
    compounds = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"integrations": [2], "multiplicities": ['s']},
            "alkyl": {"integrations": [12], "multiplicities": ['s']}
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": {"integrations": [1, 1], "multiplicities": ['s', 's']},
            "alkyl": {"integrations": [6, 3, 3], "multiplicities": ['s', 's', 's']}
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": {"integrations": [2], "multiplicities": ['s']},
            "alkyl": {"integrations": [6, 6], "multiplicities": ['s', 's']}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"integrations": [4], "multiplicities": ['s']},
            "alkyl": {"integrations": [4, 6], "multiplicities": ['q', 't']} # Quartet and Triplet
        }
    }

    # Define the options as pairs of compounds.
    options = {
        "A": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "D": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"]
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = "C"

    # Define the experimental data from the question.
    expected_aromatic = {"signals": 2, "multiplicity": "s", "ratio": [1.0, 1.0]}
    expected_alkyl = {"signals": 3, "multiplicity": "s", "ratio": [2.0, 1.0, 1.0]}

    # Iterate through each option to find the one that matches the experimental data.
    identified_correct_option = None
    for option_key, compound_names in options.items():
        c1_name, c2_name = compound_names
        c1 = compounds[c1_name]
        c2 = compounds[c2_name]

        # Combine the signals for the 1:1 mixture.
        mix_aromatic_integrations = c1["aromatic"]["integrations"] + c2["aromatic"]["integrations"]
        mix_aromatic_multiplicities = c1["aromatic"]["multiplicities"] + c2["aromatic"]["multiplicities"]
        
        mix_alkyl_integrations = c1["alkyl"]["integrations"] + c2["alkyl"]["integrations"]
        mix_alkyl_multiplicities = c1["alkyl"]["multiplicities"] + c2["alkyl"]["multiplicities"]

        # --- Check Aromatic Region ---
        # 1. Check number of signals.
        if len(mix_aromatic_integrations) != expected_aromatic["signals"]:
            continue
        
        # 2. Check if all signals are singlets.
        if not all(m == expected_aromatic["multiplicity"] for m in mix_aromatic_multiplicities):
            continue

        # 3. Check integration ratio.
        aromatic_integrations_sorted = sorted(mix_aromatic_integrations, reverse=True)
        smallest_aromatic = aromatic_integrations_sorted[-1]
        if smallest_aromatic == 0: continue
        calculated_aromatic_ratio = [i / smallest_aromatic for i in aromatic_integrations_sorted]
        if calculated_aromatic_ratio != expected_aromatic["ratio"]:
            continue

        # --- Check Alkyl Region ---
        # 1. Check number of signals.
        if len(mix_alkyl_integrations) != expected_alkyl["signals"]:
            continue

        # 2. Check if all signals are singlets.
        if not all(m == expected_alkyl["multiplicity"] for m in mix_alkyl_multiplicities):
            continue
            
        # 3. Check integration ratio.
        alkyl_integrations_sorted = sorted(mix_alkyl_integrations, reverse=True)
        smallest_alkyl = alkyl_integrations_sorted[-1]
        if smallest_alkyl == 0: continue
        calculated_alkyl_ratio = [i / smallest_alkyl for i in alkyl_integrations_sorted]
        if calculated_alkyl_ratio != expected_alkyl["ratio"]:
            continue
            
        # If all checks pass, this is the correct option.
        identified_correct_option = option_key
        break

    # --- Final Verification ---
    if identified_correct_option is None:
        return "The provided answer is incorrect. No option perfectly matches the experimental data without assuming accidental signal overlap. There might be an error in the problem statement or the analysis of individual spectra."

    if llm_answer == identified_correct_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {identified_correct_option}. "
                f"The mixture in option {identified_correct_option} ({' and '.join(options[identified_correct_option])}) "
                f"is the only one that produces {expected_aromatic['signals']} aromatic singlets in a 1:1 ratio and "
                f"{expected_alkyl['signals']} alkyl singlets in a 2:1:1 ratio.")

# Run the check
result = check_nmr_answer()
print(result)