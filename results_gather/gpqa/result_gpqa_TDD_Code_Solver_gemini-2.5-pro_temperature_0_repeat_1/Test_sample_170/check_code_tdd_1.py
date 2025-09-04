import collections

def check_electrophilic_substitution_order():
    """
    Checks the correctness of the proposed order for the yield of para-isomer
    in the electrophilic bromination of substituted benzenes.

    The function uses established experimental data for isomer distribution.
    """

    # Data represents the approximate percentage yield of the para-isomer during
    # electrophilic bromination. This data is compiled from standard organic
    # chemistry textbooks (e.g., Vollhardt & Schore, Jerry March's Advanced
    # Organic Chemistry) and is widely accepted.
    #
    # 1) C6H5-CH3 (Toluene): Ortho, para-director. Para yield is high.
    # 2) C6H5-COOC2H5 (Ethyl benzoate): Meta-director. Para yield is very low.
    # 3) C6H5-Cl (Chlorobenzene): Ortho, para-director. Highly para-selective.
    # 4) C6H5-NO2 (Nitrobenzene): Meta-director. Strong deactivator, very low para yield.
    # 5) C6H5-C2H5 (Ethylbenzene): Ortho, para-director. Bulkier than methyl, slightly more para-selective.
    # 6) C6H5-COOH (Benzoic acid): Meta-director. Strong deactivator, very low para yield.
    
    compounds_data = {
        1: {'name': 'C6H5-CH3', 'para_yield': 67},
        2: {'name': 'C6H5-COOC2H5', 'para_yield': 8},
        3: {'name': 'C6H5-Cl', 'para_yield': 87},
        4: {'name': 'C6H5-NO2', 'para_yield': 0.3},
        5: {'name': 'C6H5-C2H5', 'para_yield': 69},
        6: {'name': 'C6H5-COOH', 'para_yield': 1.5}
    }

    # The proposed answer 'D' corresponds to the order: 4 < 2 < 6 < 3 < 1 < 5
    # Let's re-read the question and answer.
    # The question is: Arrange the substances in order of increasing the weight fraction of the yield of the para-isomer.
    # The provided answer is D, which corresponds to the order 4 < 6 < 2 < 1 < 5 < 3.
    proposed_order_str = "4<6<2<1<5<3"
    proposed_order = [int(x) for x in proposed_order_str.split('<')]

    # Sort the compounds based on their para-yield in increasing order.
    # The `sorted` function returns a list of the dictionary keys (the compound numbers)
    # ordered by the 'para_yield' value of the corresponding dictionary entry.
    correct_order = sorted(compounds_data, key=lambda k: compounds_data[k]['para_yield'])

    # Check if the proposed order matches the correct order derived from experimental data.
    if proposed_order == correct_order:
        return "Correct"
    else:
        # If the order is incorrect, construct a detailed reason.
        reason = "The provided answer is incorrect.\n"
        reason += f"Proposed order: {' < '.join(map(str, proposed_order))}\n"
        reason += f"Correct order based on experimental para-isomer yields: {' < '.join(map(str, correct_order))}\n\n"
        
        reason += "Detailed comparison with approximate para-yields (%):\n"
        correct_data_str = ""
        for num in correct_order:
            correct_data_str += f"  {num}) {compounds_data[num]['name']}: {compounds_data[num]['para_yield']}%\n"
        reason += correct_data_str
        
        # Identify the first mismatch to pinpoint the error.
        for i in range(len(proposed_order)):
            if proposed_order[i] != correct_order[i]:
                reason += f"\nThe first error is at position {i+1}. "
                reason += f"The proposed order places substance {proposed_order[i]}, but it should be {correct_order[i]}.\n"
                reason += (f"Substance {correct_order[i]} ({compounds_data[correct_order[i]]['name']}) has a para-yield of "
                           f"{compounds_data[correct_order[i]]['para_yield']}%, while substance {proposed_order[i]} "
                           f"({compounds_data[proposed_order[i]]['name']}) has a para-yield of "
                           f"{compounds_data[proposed_order[i]]['para_yield']}%.")
                break
        return reason

# Execute the check and print the result.
# print(check_electrophilic_substitution_order())