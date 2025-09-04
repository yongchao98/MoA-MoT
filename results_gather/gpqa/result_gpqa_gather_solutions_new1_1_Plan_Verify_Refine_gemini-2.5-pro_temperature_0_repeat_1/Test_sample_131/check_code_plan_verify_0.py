import math

def check_nmr_answer():
    """
    Checks the correctness of the answer for the given 1H NMR spectroscopy problem.
    It models the predicted spectra for each compound and checks which 1:1 mixture
    matches the experimental data.
    """

    # Step 1: Define the predicted 1H NMR data for each compound based on its structure.
    # Integrals represent the number of protons for each signal.
    compounds_data = {
        "1,2,4,5-tetramethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrals": [2]},
            "aliphatic": {"multiplicity": ["singlet"], "integrals": [12]}
        },
        "1,2,3,5-tetramethylbenzene": {
            # Has two non-equivalent aromatic protons (H4, H6) -> 2 singlets
            "aromatic": {"multiplicity": ["singlet", "singlet"], "integrals": [1, 1]},
            # Has three sets of methyl groups (C1/C3, C2, C5) -> 3 singlets
            "aliphatic": {"multiplicity": ["singlet", "singlet", "singlet"], "integrals": [6, 3, 3]}
        },
        "1,2,3,4-tetramethylbenzene": {
            # Aromatic protons at C5/C6 form an AA' system, often observed as a singlet.
            "aromatic": {"multiplicity": ["singlet"], "integrals": [2]},
            # Two sets of methyl groups (C1/C4, C2/C3) -> 2 singlets
            "aliphatic": {"multiplicity": ["singlet", "singlet"], "integrals": [6, 6]}
        },
        "1,4-diethylbenzene": {
            "aromatic": {"multiplicity": ["singlet"], "integrals": [4]},
            # Ethyl groups give a quartet and a triplet.
            "aliphatic": {"multiplicity": ["quartet", "triplet"], "integrals": [4, 6]}
        }
    }

    # Step 2: Define the options from the question.
    options = {
        "A": ["1,2,4,5-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "B": ["1,2,4,5-tetramethylbenzene", "1,2,3,4-tetramethylbenzene"],
        "C": ["1,2,3,4-tetramethylbenzene", "1,2,3,5-tetramethylbenzene"],
        "D": ["1,2,3,5-tetramethylbenzene", "1,4-diethylbenzene"]
    }

    # The final answer provided by the LLM to be checked.
    llm_answer = "B"

    # Step 3: Define a function to check if a given option matches the experimental data.
    def check_option(option_letter):
        c1_name, c2_name = options[option_letter]
        c1 = compounds_data[c1_name]
        c2 = compounds_data[c2_name]

        # Constraint: All signals must be singlets.
        for m in c1["aliphatic"]["multiplicity"] + c2["aliphatic"]["multiplicity"]:
            if m != "singlet":
                return f"Option {option_letter} is incorrect. It contains '{c1_name if m in c1['aliphatic']['multiplicity'] else c2_name}', which has a non-singlet ('{m}') in the aliphatic region."
        
        # Constraint: Aromatic region must have 2 singlets in a 1:1 ratio.
        aromatic_signals_count = len(c1["aromatic"]["integrals"]) + len(c2["aromatic"]["integrals"])
        if aromatic_signals_count != 2:
            return f"Option {option_letter} is incorrect. The mixture is predicted to have {aromatic_signals_count} aromatic signals, but the question specifies 2."
        
        aromatic_integrals = c1["aromatic"]["integrals"] + c2["aromatic"]["integrals"]
        if aromatic_integrals[0] != aromatic_integrals[1]:
            return f"Option {option_letter} is incorrect. The aromatic signals have integrals {aromatic_integrals}, which is not a 1:1 ratio."

        # Constraint: Aliphatic region must have 3 singlets in a 2:1:1 ratio.
        aliphatic_signals_count = len(c1["aliphatic"]["integrals"]) + len(c2["aliphatic"]["integrals"])
        if aliphatic_signals_count != 3:
            return f"Option {option_letter} is incorrect. The mixture is predicted to have {aliphatic_signals_count} aliphatic signals, but the question specifies 3."

        aliphatic_integrals = sorted(c1["aliphatic"]["integrals"] + c2["aliphatic"]["integrals"], reverse=True)
        # For a 2:1:1 ratio, the largest integral must be twice the second, and the second and third must be equal.
        if not (aliphatic_integrals[0] == 2 * aliphatic_integrals[1] and aliphatic_integrals[1] == aliphatic_integrals[2]):
            return f"Option {option_letter} is incorrect. The aliphatic signals have integrals {aliphatic_integrals}, which is not a 2:1:1 ratio."
            
        return "Correct"

    # Step 4: Evaluate the LLM's answer and check for other possibilities.
    results = {opt: check_option(opt) for opt in options}
    
    if results.get(llm_answer) == "Correct":
        correct_options = [opt for opt, res in results.items() if res == "Correct"]
        if len(correct_options) == 1:
            return "Correct"
        else:
            # This case would occur if multiple options satisfied the criteria.
            return f"Ambiguous: The provided answer '{llm_answer}' is correct, but other options { [opt for opt in correct_options if opt != llm_answer] } also satisfy the conditions."
    else:
        # The provided answer is incorrect. State the reason.
        reason = results.get(llm_answer, "is not a valid option.")
        return f"Incorrect. The provided answer '{llm_answer}' is wrong. Reason: {reason}"

# Execute the check and print the result.
result = check_nmr_answer()
print(result)