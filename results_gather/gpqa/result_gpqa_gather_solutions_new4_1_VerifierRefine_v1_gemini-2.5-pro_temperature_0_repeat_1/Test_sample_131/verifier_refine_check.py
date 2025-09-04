import collections

def check_correctness():
    """
    Checks the correctness of the LLM's answer by simulating the 1H NMR spectrum
    for each possible mixture and comparing it to the experimental data.
    """

    # Ground truth NMR data for each compound based on chemical principles.
    # Signal format: list of tuples (integration, multiplicity)
    compounds_db = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": [(2, "singlet")],
            "aliphatic": [(12, "singlet")]
        },
        "1,2,3,5-tetramethylbenzene": {
            "aromatic": [(1, "singlet"), (1, "singlet")],
            "aliphatic": [(6, "singlet"), (3, "singlet"), (3, "singlet")]
        },
        "1,2,3,4-tetramethylbenzene": {
            "aromatic": [(2, "singlet")],
            "aliphatic": [(6, "singlet"), (6, "singlet")]
        },
        "1,4-diethylbenzene": {
            "aromatic": [(4, "singlet")],
            "aliphatic": [(4, "quartet"), (6, "triplet")]
        }
    }

    # Define the options as given in the question
    options = {
        "A": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "C": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"]
    }

    # Define the target spectrum from the problem description
    target_spectrum = {
        "aromatic": {
            "count": 2,
            "multiplicity": "singlet",
            "ratio": sorted([1, 1])
        },
        "aliphatic": {
            "count": 3,
            "multiplicity": "singlet",
            "ratio": sorted([2, 1, 1])
        }
    }

    # The final answer provided by the LLM to be checked
    llm_answer = "B"

    # --- Verification Logic ---
    
    def get_ratio(integrations):
        """Normalizes a list of integrations into a ratio."""
        if not integrations:
            return []
        integrations.sort()
        smallest = integrations[0]
        # Handle potential division by zero, though unlikely here
        if smallest == 0:
            return [0] * len(integrations)
        ratio = [int(round(i / smallest)) for i in integrations]
        return ratio

    correct_option_found = None

    for option_letter, compound_names in options.items():
        comp1_name, comp2_name = compound_names
        comp1_data = compounds_db[comp1_name]
        comp2_data = compounds_db[comp2_name]

        # --- Combine signals for the 1:1 mixture ---
        
        # Aromatic signals
        aromatic_signals = comp1_data["aromatic"] + comp2_data["aromatic"]
        aromatic_integrations = [integ for integ, _ in aromatic_signals]
        aromatic_multiplicities = [mult for _, mult in aromatic_signals]

        # Aliphatic signals
        aliphatic_signals = comp1_data["aliphatic"] + comp2_data["aliphatic"]
        aliphatic_integrations = [integ for integ, _ in aliphatic_signals]
        aliphatic_multiplicities = [mult for _, mult in aliphatic_signals]

        # --- Check against experimental data ---
        
        # Check aromatic signal count
        if len(aromatic_signals) != target_spectrum["aromatic"]["count"]:
            continue
        
        # Check aromatic multiplicity (all must be singlets)
        if not all(m == target_spectrum["aromatic"]["multiplicity"] for m in aromatic_multiplicities):
            continue
            
        # Check aromatic ratio
        if get_ratio(aromatic_integrations) != target_spectrum["aromatic"]["ratio"]:
            continue

        # Check aliphatic signal count
        if len(aliphatic_signals) != target_spectrum["aliphatic"]["count"]:
            continue

        # Check aliphatic multiplicity (all must be singlets)
        if not all(m == target_spectrum["aliphatic"]["multiplicity"] for m in aliphatic_multiplicities):
            continue

        # Check aliphatic ratio
        if get_ratio(aliphatic_integrations) != target_spectrum["aliphatic"]["ratio"]:
            continue
            
        # If all conditions are met, this is the correct option
        correct_option_found = option_letter
        break

    # --- Final Verdict ---
    if correct_option_found == llm_answer:
        return "Correct"
    else:
        if correct_option_found is None:
            return f"Incorrect. The provided answer is {llm_answer}, but no option matches the experimental data based on the analysis. There might be an error in the problem statement or the options."
        else:
            return (f"Incorrect. The provided answer is {llm_answer}, but the correct answer is {correct_option_found}. "
                    f"The mixture in option {correct_option_found} ({', '.join(options[correct_option_found])}) is the only one that matches all spectral data. "
                    f"Aromatic check: {options[correct_option_found][0]} gives a 2H singlet and {options[correct_option_found][1]} gives a 2H singlet, resulting in two 1:1 signals. "
                    f"Aliphatic check: {options[correct_option_found][0]} gives a 12H singlet and {options[correct_option_found][1]} gives two 6H singlets, resulting in three signals with a 12:6:6 or 2:1:1 ratio.")

# Execute the check and print the result
result = check_correctness()
print(result)