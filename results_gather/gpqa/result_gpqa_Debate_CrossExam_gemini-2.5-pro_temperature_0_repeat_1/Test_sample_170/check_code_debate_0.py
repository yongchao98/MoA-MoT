def check_electrophilic_substitution_order():
    """
    Checks the correctness of the given order for the yield of para-isomer in electrophilic bromination.

    The function uses established chemical principles to rank the substances:
    1.  Meta-directors yield less para-isomer than ortho/para-directors.
    2.  Within meta-directors, stronger deactivation leads to lower para-yield.
        Order of deactivating strength: -NO2 > -COOH > -COOC2H5
        Therefore, order of increasing para-yield: -NO2 < -COOH < -COOC2H5 (4 < 6 < 2)
    3.  Within ortho/para-directors, steric and electronic effects determine the para/ortho ratio.
        - Sterics: The bulkier -C2H5 group favors para more than -CH3. So, Toluene < Ethylbenzene (1 < 5).
        - Electronics: The -Cl group, due to its strong inductive effect at the ortho position,
          is highly para-directing, more so than alkyl groups. So, Ethylbenzene < Chlorobenzene (5 < 3).
        - Therefore, order of increasing para-yield: -CH3 < -C2H5 < -Cl (1 < 5 < 3).
    4.  Combining these gives the final correct order.
    """

    # Define the substances and their properties relevant to para-isomer yield.
    # 'para_yield_rank' is a numerical representation of the expected yield, where a higher number means a higher yield.
    substances = {
        1: {'name': 'Toluene (C6H5-CH3)', 'rank': 4},
        2: {'name': 'Ethyl benzoate (C6H5-COOC2H5)', 'rank': 3},
        3: {'name': 'Chlorobenzene (C6H5-Cl)', 'rank': 6},
        4: {'name': 'Nitrobenzene (C6H5-NO2)', 'rank': 1},
        5: {'name': 'Ethylbenzene (C6H5-C2H5)', 'rank': 5},
        6: {'name': 'Benzoic acid (C6H5-COOH)', 'rank': 2}
    }

    # The answer provided by the LLM is B, which corresponds to the order: 4 < 6 < 2 < 1 < 5 < 3
    llm_answer_string = "B"
    llm_answer_map = {
        "A": [4, 2, 6, 3, 1, 5],
        "B": [4, 6, 2, 1, 5, 3],
        "C": [3, 5, 1, 6, 2, 4],
        "D": [6, 2, 4, 5, 1, 3]
    }

    if llm_answer_string not in llm_answer_map:
        return f"The provided answer '{llm_answer_string}' is not a valid option."

    llm_order = llm_answer_map[llm_answer_string]

    # Determine the correct order by sorting the substances based on their rank.
    # The `sorted` function takes the dictionary items and sorts them based on the 'rank' value.
    sorted_substances = sorted(substances.items(), key=lambda item: item[1]['rank'])
    
    # Extract the substance numbers to get the correct sequence.
    correct_order = [item[0] for item in sorted_substances]

    # Compare the LLM's answer with the correct order.
    if llm_order == correct_order:
        return "Correct"
    else:
        # If the order is incorrect, provide a detailed reason.
        reason = f"The answer {llm_answer_string} is incorrect.\n"
        reason += f"Proposed order: {llm_order}\n"
        reason += f"Correct order based on chemical principles: {correct_order}\n\n"
        reason += "Detailed explanation for the correct order (4 < 6 < 2 < 1 < 5 < 3):\n"
        reason += "1. Meta-directors (-NO2, -COOH, -COOC2H5) yield the least para-product. Their increasing para-yield order is based on decreasing deactivating strength: NO2 (4) < COOH (6) < COOC2H5 (2).\n"
        reason += "2. Ortho/para-directors (-CH3, -C2H5, -Cl) yield the most para-product. Their increasing para-yield order is:\n"
        reason += "   - Toluene (1) has a lower para/ortho ratio than Ethylbenzene (5) due to lower steric hindrance.\n"
        reason += "   - Ethylbenzene (5) has a lower para/ortho ratio than Chlorobenzene (3) because the strong inductive effect of Cl disfavors ortho attack more significantly, making para-substitution highly favorable.\n"
        reason += "   - This gives the order: 1 < 5 < 3.\n"
        reason += "3. Combining both groups results in the final order: 4 < 6 < 2 < 1 < 5 < 3."
        return reason

# Execute the check and print the result.
result = check_electrophilic_substitution_order()
print(result)