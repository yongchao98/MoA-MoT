import collections

def check_correctness():
    """
    Checks the correctness of the answer to the NMR spectroscopy question.

    This function simulates the analysis of 1H NMR spectra for mixtures of
    aromatic compounds to verify the given answer.
    """

    # Predicted 1H NMR data for each individual compound.
    # Format: { 'aromatic': [(multiplicity, integration), ...], 'alkyl': [...] }
    spectra = {
        '1,2,4,5-tetramethylbenzene': {
            'aromatic': [('singlet', 2)],
            'alkyl': [('singlet', 12)]
        },
        '1,2,3,5-tetramethylbenzene': {
            'aromatic': [('singlet', 1), ('singlet', 1)],
            'alkyl': [('singlet', 6), ('singlet', 3), ('singlet', 3)]
        },
        '1,2,3,4-tetramethylbenzene': {
            'aromatic': [('singlet', 2)],
            'alkyl': [('singlet', 6), ('singlet', 6)]
        },
        '1,4-diethylbenzene': {
            'aromatic': [('singlet', 4)],
            'alkyl': [('quartet', 4), ('triplet', 6)]
        }
    }

    # Define the options from the question
    options = {
        'A': ['1,2,3,4-tetramethylbenzene', '1,2,3,5-tetramethylbenzene'],
        'B': ['1,2,4,5-tetramethylbenzene', '1,2,3,4-tetramethylbenzene'],
        'C': ['1,2,4,5-tetramethylbenzene', '1,2,3,5-tetramethylbenzene'],
        'D': ['1,2,3,5-tetramethylbenzene', '1,4-diethylbenzene']
    }

    # The final answer provided by the LLM
    llm_answer = 'B'

    # --- Helper function to check a single mixture ---
    def check_mixture(compound1_name, compound2_name):
        c1_spec = spectra[compound1_name]
        c2_spec = spectra[compound2_name]

        # Combine the signals from both compounds for a 1:1 mixture
        combined_aromatic = c1_spec['aromatic'] + c2_spec['aromatic']
        combined_alkyl = c1_spec['alkyl'] + c2_spec['alkyl']

        # --- Check Aromatic Region Constraints ---
        # Constraint: Two signals
        if len(combined_aromatic) != 2:
            return False, f"Fails aromatic constraint: Expected 2 signals, but found {len(combined_aromatic)}."
        # Constraint: All singlets
        if not all(sig[0] == 'singlet' for sig in combined_aromatic):
            return False, "Fails aromatic constraint: Not all signals are singlets."
        # Constraint: 1:1 integration ratio
        aromatic_integrations = [sig[1] for sig in combined_aromatic]
        if aromatic_integrations[0] != aromatic_integrations[1]:
            return False, f"Fails aromatic constraint: Integration ratio is not 1:1 (found {aromatic_integrations})."

        # --- Check Alkyl Region Constraints ---
        # Constraint: Three signals
        if len(combined_alkyl) != 3:
            return False, f"Fails alkyl constraint: Expected 3 signals, but found {len(combined_alkyl)}."
        # Constraint: All singlets
        if not all(sig[0] == 'singlet' for sig in combined_alkyl):
            return False, "Fails alkyl constraint: Not all signals are singlets."
        # Constraint: 2:1:1 integration ratio
        alkyl_integrations = sorted([sig[1] for sig in combined_alkyl], reverse=True)
        # Check the ratio by dividing by the smallest integration value
        if not (alkyl_integrations[0] / alkyl_integrations[2] == 2.0 and \
                alkyl_integrations[1] / alkyl_integrations[2] == 1.0 and \
                alkyl_integrations[2] / alkyl_integrations[2] == 1.0):
            return False, f"Fails alkyl constraint: Integration ratio is not 2:1:1 (found {alkyl_integrations})."

        return True, "Matches all constraints."

    # --- Evaluate all options ---
    matching_options = []
    for option_letter, compounds in options.items():
        is_match, reason = check_mixture(compounds[0], compounds[1])
        if is_match:
            matching_options.append(option_letter)

    # --- Determine final correctness ---
    if llm_answer in matching_options and len(matching_options) == 1:
        return "Correct"
    elif llm_answer not in matching_options:
        # Find the reason why the LLM's answer failed
        _, reason_for_failure = check_mixture(options[llm_answer][0], options[llm_answer][1])
        return f"Incorrect. The proposed answer {llm_answer} is wrong. {reason_for_failure}"
    elif len(matching_options) > 1:
        return f"Incorrect. The proposed answer {llm_answer} is a valid solution, but other options also match: {matching_options}. The question may be ambiguous."
    else: # This case should not be reached if there is a correct answer
        return "Incorrect. No option correctly matches the given NMR data."

# Run the check and print the result
result = check_correctness()
print(result)