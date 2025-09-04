def check_electrophilic_substitution_order():
    """
    Checks the correctness of the given answer for arranging substances by increasing para-isomer yield in bromination.
    """
    
    # The provided answer from the other LLM to be checked.
    llm_answer = "A"
    
    # Map of options to their corresponding numerical sequences for the substances 1-6.
    options = {
        "A": [4, 6, 2, 1, 5, 3],
        "B": [3, 5, 1, 6, 2, 4],
        "C": [6, 2, 4, 5, 1, 3],
        "D": [4, 2, 6, 3, 1, 5]
    }

    if llm_answer not in options:
        return f"The provided answer '{llm_answer}' is not a valid option (A, B, C, or D)."

    llm_sequence = options[llm_answer]

    # This data structure encodes the chemical principles. The 'para_yield' is a representative
    # value based on experimental data and common textbook rules, used here for sorting.
    # 1: Toluene, 2: Ethyl benzoate, 3: Chlorobenzene, 4: Nitrobenzene, 5: Ethylbenzene, 6: Benzoic acid
    substance_data = {
        1: {'name': 'Toluene (C6H5-CH3)', 'para_yield': 67, 'type': 'o,p-director'},
        2: {'name': 'Ethyl benzoate (C6H5-COOC2H5)', 'para_yield': 6, 'type': 'meta-director'},
        3: {'name': 'Chlorobenzene (C6H5-Cl)', 'para_yield': 87, 'type': 'o,p-director'},
        4: {'name': 'Nitrobenzene (C6H5-NO2)', 'para_yield': 0.3, 'type': 'meta-director'},
        5: {'name': 'Ethylbenzene (C6H5-C2H5)', 'para_yield': 75, 'type': 'o,p-director'},
        6: {'name': 'Benzoic acid (C6H5-COOH)', 'para_yield': 1.5, 'type': 'meta-director'}
    }

    # Determine the correct order by sorting the substances based on their para_yield.
    correct_sequence = sorted(substance_data.keys(), key=lambda k: substance_data[k]['para_yield'])

    # Check if the LLM's sequence matches the derived correct sequence.
    if llm_sequence == correct_sequence:
        return "Correct"
    else:
        correct_order_str = '<'.join(map(str, correct_sequence))
        llm_order_str = '<'.join(map(str, llm_sequence))
        
        reason = f"The answer '{llm_answer}' corresponding to the sequence {llm_order_str} is incorrect.\n"
        reason += f"The correct order of increasing para-isomer yield is {correct_order_str}.\n\n"
        reason += "Justification:\n"
        reason += "1. Meta-directing groups (4, 6, 2) produce the least para-isomer. Their order based on increasing yield is 4(-NO2) < 6(-COOH) < 2(-COOC2H5).\n"
        reason += "2. Ortho,para-directing groups (1, 5, 3) produce the most para-isomer. Their order based on increasing yield is 1(-CH3) < 5(-C2H5) < 3(-Cl), due to steric hindrance and electronic effects.\n"
        reason += f"3. Combining these gives the final correct order: {correct_order_str}."

        return reason

# Execute the check and print the result.
print(check_electrophilic_substitution_order())