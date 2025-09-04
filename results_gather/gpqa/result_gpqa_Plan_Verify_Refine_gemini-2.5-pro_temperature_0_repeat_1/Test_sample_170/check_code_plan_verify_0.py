def check_answer():
    """
    Checks the correctness of the given answer for arranging substances by increasing para-isomer yield in bromination.
    """
    # The answer to be checked. Option B corresponds to the order [4, 6, 2, 1, 5, 3]
    given_answer_order = [4, 6, 2, 1, 5, 3]

    # Define the properties of each substance based on chemical principles.
    # 'type': 'meta' or 'op' (ortho,para)
    # 'rank': A numerical rank for sorting. Lower rank means lower para-isomer yield.
    # The ranks are assigned based on established chemical principles:
    # 1. Meta-directors have lower para-yield than o,p-directors.
    # 2. For meta-directors, yield increases as deactivating strength decreases: NO2 < COOH < COOC2H5.
    # 3. For o,p-directors, yield increases with steric hindrance (CH3 < C2H5) and is highest for Cl.
    substances = {
        1: {'name': 'Toluene (C6H5-CH3)', 'type': 'op', 'rank': 4},
        2: {'name': 'Ethyl benzoate (C6H5-COOC2H5)', 'type': 'meta', 'rank': 3},
        3: {'name': 'Chlorobenzene (C6H5-Cl)', 'type': 'op', 'rank': 6},
        4: {'name': 'Nitrobenzene (C6H5-NO2)', 'type': 'meta', 'rank': 1},
        5: {'name': 'Ethylbenzene (C6H5-C2H5)', 'type': 'op', 'rank': 5},
        6: {'name': 'Benzoic acid (C6H5-COOH)', 'type': 'meta', 'rank': 2}
    }

    # Sort the substance keys based on their assigned rank
    # The key for sorting is a lambda function that retrieves the 'rank' for each substance number.
    correct_order = sorted(substances.keys(), key=lambda x: substances[x]['rank'])

    # Check if the generated correct order matches the given answer's order
    if correct_order == given_answer_order:
        return "Correct"
    else:
        # If they don't match, construct a detailed error message.
        error_message = "The answer is incorrect. The calculated correct order is {}.\n".format(correct_order)
        error_message += "The provided answer's order is {}.\n".format(given_answer_order)
        
        # Find the first point of deviation
        for i in range(len(correct_order)):
            if correct_order[i] != given_answer_order[i]:
                correct_item = substances[correct_order[i]]['name']
                given_item = substances[given_answer_order[i]]['name']
                error_message += "The first discrepancy is at position {}: ".format(i+1)
                error_message += "The correct substance should be {} but the answer has {}.\n".format(correct_item, given_item)
                break
        
        error_message += "\nReasoning for the correct order (4 < 6 < 2 < 1 < 5 < 3):\n"
        error_message += "1. Meta-directors (-NO2, -COOH, -COOC2H5) have the lowest para-yield.\n"
        error_message += "   - Their order is based on deactivating strength (para-yield: NO2 < COOH < COOC2H5), so the sequence starts with 4 < 6 < 2.\n"
        error_message += "2. Ortho,para-directors (-CH3, -C2H5, -Cl) have higher para-yield.\n"
        error_message += "   - Steric hindrance from the bulkier -C2H5 group favors the para position more than -CH3 does. So, para-yield: Toluene(1) < Ethylbenzene(5).\n"
        error_message += "   - Chlorobenzene(3) gives a very high proportion of the para-isomer, the highest among these three.\n"
        error_message += "   - This gives the sequence 1 < 5 < 3 for this group."

        return error_message

# Execute the check and print the result
result = check_answer()
print(result)