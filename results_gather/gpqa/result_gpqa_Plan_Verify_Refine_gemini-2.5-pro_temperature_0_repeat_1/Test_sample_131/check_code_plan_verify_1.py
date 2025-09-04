import math
from functools import reduce

def check_correctness_of_nmr_answer():
    """
    This function checks the correctness of the provided LLM answer by systematically
    analyzing the 1H NMR spectral data for each possible mixture of compounds.
    """

    # Step 1: Define the predicted 1H NMR data for each individual compound.
    # The data includes the number, multiplicity, and integration of signals in both
    # the aromatic and aliphatic regions.
    # Structure: { "aromatic_signals": [(multiplicity, integration), ...], "aliphatic_signals": [...] }
    compounds_data = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic_signals": [("singlet", 2)],  # One signal for 2 equivalent aromatic H's
            "aliphatic_signals": [("singlet", 12)] # One signal for 4 equivalent methyl groups
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic_signals": [("singlet", 1), ("singlet", 1)], # Two signals for 2 non-equivalent aromatic H's
            "aliphatic_signals": [("singlet", 6), ("singlet", 3), ("singlet", 3)] # Three signals for 3 sets of non-equivalent methyl groups
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic_signals": [("singlet", 2)], # One signal for 2 equivalent aromatic H's
            "aliphatic_signals": [("singlet", 6), ("singlet", 6)] # Two signals for 2 sets of equivalent methyl groups
        },
        "1,4-diethylbenzene": {
            "aromatic_signals": [("singlet", 4)], # One signal for 4 equivalent aromatic H's
            "aliphatic_signals": [("quartet", 4), ("triplet", 6)] # A quartet (CH2) and a triplet (CH3)
        }
    }

    # Step 2: Define the experimental data given in the question.
    # Ratios are sorted for consistent comparison.
    experimental_data = {
        "aromatic": {
            "num_signals": 2,
            "multiplicities": {"singlet"},
            "ratio": sorted([1, 1], reverse=True)
        },
        "aliphatic": {
            "num_signals": 3,
            "multiplicities": {"singlet"},
            "ratio": sorted([2, 1, 1], reverse=True)
        }
    }

    # Step 3: Define the multiple-choice options and the LLM's answer.
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "B": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "D": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"]
    }
    llm_answer = "C"

    # Helper function to simplify a list of integrations into their smallest integer ratio.
    def get_simplified_ratio(integrations):
        if not integrations:
            return []
        if len(integrations) == 1:
            return [1]
        common_divisor = reduce(math.gcd, integrations)
        ratio = [i // common_divisor for i in integrations]
        return sorted(ratio, reverse=True)

    # Step 4: Analyze each option against the experimental data.
    results = {}
    for option, compound_names in options.items():
        comp1_name, comp2_name = compound_names
        comp1 = compounds_data[comp1_name]
        comp2 = compounds_data[comp2_name]

        # For a 1:1 mixture, we simply combine the signals from both compounds.
        combined_aromatic = comp1["aromatic_signals"] + comp2["aromatic_signals"]
        combined_aliphatic = comp1["aliphatic_signals"] + comp2["aliphatic_signals"]

        # --- Check Aromatic Region ---
        if len(combined_aromatic) != experimental_data["aromatic"]["num_signals"]:
            results[option] = f"Fails on aromatic signal count. Expected {experimental_data['aromatic']['num_signals']}, but mixture gives {len(combined_aromatic)}."
            continue
        
        aromatic_multiplicities = {sig[0] for sig in combined_aromatic}
        if aromatic_multiplicities != experimental_data["aromatic"]["multiplicities"]:
            results[option] = f"Fails on aromatic multiplicity. Expected only {experimental_data['aromatic']['multiplicities']}, but found {aromatic_multiplicities}."
            continue

        aromatic_integrations = [sig[1] for sig in combined_aromatic]
        calculated_aromatic_ratio = get_simplified_ratio(aromatic_integrations)
        if calculated_aromatic_ratio != experimental_data["aromatic"]["ratio"]:
            results[option] = f"Fails on aromatic ratio. Expected {experimental_data['aromatic']['ratio']}, but calculated {calculated_aromatic_ratio}."
            continue

        # --- Check Aliphatic Region ---
        if len(combined_aliphatic) != experimental_data["aliphatic"]["num_signals"]:
            results[option] = f"Fails on aliphatic signal count. Expected {experimental_data['aliphatic']['num_signals']}, but mixture gives {len(combined_aliphatic)}."
            continue

        aliphatic_multiplicities = {sig[0] for sig in combined_aliphatic}
        if aliphatic_multiplicities != experimental_data["aliphatic"]["multiplicities"]:
            results[option] = f"Fails on aliphatic multiplicity. Expected only {experimental_data['aliphatic']['multiplicities']}, but found {aliphatic_multiplicities}."
            continue

        aliphatic_integrations = [sig[1] for sig in combined_aliphatic]
        calculated_aliphatic_ratio = get_simplified_ratio(aliphatic_integrations)
        if calculated_aliphatic_ratio != experimental_data["aliphatic"]["ratio"]:
            results[option] = f"Fails on aliphatic ratio. Expected {experimental_data['aliphatic']['ratio']}, but calculated {calculated_aliphatic_ratio} from integrations {aliphatic_integrations}."
            continue
            
        # If all checks pass for this option
        results[option] = "Pass"

    # Step 5: Evaluate the LLM's answer based on the analysis.
    passing_options = [opt for opt, res in results.items() if res == "Pass"]

    if llm_answer in passing_options:
        if len(passing_options) == 1:
            # The LLM's answer is the only one that passed all checks.
            return "Correct"
        else:
            # This case is unlikely for this problem but included for robustness.
            other_correct = [opt for opt in passing_options if opt != llm_answer]
            return f"The answer {llm_answer} is correct, but the analysis is incomplete as option(s) {other_correct} also satisfy the conditions."
    else:
        # The LLM's answer is incorrect.
        reason = results.get(llm_answer, "The option was not evaluated or does not exist.")
        if not passing_options:
            return f"The answer {llm_answer} is incorrect. Reason: {reason}. Furthermore, no other option satisfies all the conditions."
        else:
            return f"The answer {llm_answer} is incorrect. Reason: {reason}. The correct option should be {passing_options[0]}, which satisfies all conditions."

# Execute the check and print the result.
result = check_correctness_of_nmr_answer()
print(result)