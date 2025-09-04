def check_electrophilic_substitution_order():
    """
    Checks the correctness of the given answer for arranging substances by increasing para-isomer yield in bromination.

    The function verifies the answer based on established chemical principles of electrophilic aromatic substitution:
    1. Meta-directors yield negligible para-isomer compared to ortho,para-directors.
    2. Within meta-directors, yield is inversely proportional to deactivating strength (-NO2 < -COOH < -COOC2H5).
    3. Within ortho,para-directors, yield depends on para-selectivity, which increases with steric hindrance and due to electronic effects of halogens (-CH3 < -C2H5 < -Cl).
    """
    # The provided answer from the LLM is D, which corresponds to the sequence 4<6<2<1<5<3
    llm_answer_option = "D"
    llm_answer_sequence_str = "4<6<2<1<5<3"

    # Step 1: Define the correct order based on chemical principles.
    # Meta-directors ranked by increasing para-yield (inversely by deactivating strength):
    # -NO2 (4) < -COOH (6) < -COOC2H5 (2)
    meta_directors_correct_order = [4, 6, 2]

    # Ortho,para-directors ranked by increasing para-yield (by para-selectivity):
    # -CH3 (1) < -C2H5 (5) < -Cl (3)
    op_directors_correct_order = [1, 5, 3]

    # The final correct order is the sequence of meta-directors followed by ortho,para-directors.
    expected_order = meta_directors_correct_order + op_directors_correct_order

    # Step 2: Parse the sequence from the LLM's answer.
    try:
        llm_order = [int(x) for x in llm_answer_sequence_str.split('<')]
    except (ValueError, AttributeError):
        return f"The answer format is incorrect. Could not parse the sequence '{llm_answer_sequence_str}'."

    # Step 3: Check if the parsed sequence matches the expected correct order.
    if llm_order == expected_order:
        return "Correct"
    else:
        # Provide a reason for the incorrectness.
        # Check 1: Grouping (all meta before all o,p)
        meta_indices = [llm_order.index(i) for i in meta_directors_correct_order]
        op_indices = [llm_order.index(i) for i in op_directors_correct_order]
        if max(meta_indices) > min(op_indices):
            return (f"Incorrect order. A meta-directing substance is placed after an ortho,para-directing one. "
                    f"All meta-directors ({meta_directors_correct_order}) must have lower para-yields than all "
                    f"ortho,para-directors ({op_directors_correct_order}). The proposed order was {llm_order}.")

        # Check 2: Order within meta-directors
        llm_meta_order = [i for i in llm_order if i in meta_directors_correct_order]
        if llm_meta_order != meta_directors_correct_order:
            return (f"Incorrect order. The relative ranking of meta-directors is wrong. "
                    f"The correct order based on increasing para-yield is {meta_directors_correct_order} "
                    f"(-NO2 < -COOH < -COOC2H5). The proposed order for these was {llm_meta_order}.")

        # Check 3: Order within o,p-directors
        llm_op_order = [i for i in llm_order if i in op_directors_correct_order]
        if llm_op_order != op_directors_correct_order:
            return (f"Incorrect order. The relative ranking of ortho,para-directors is wrong. "
                    f"The correct order based on increasing para-yield is {op_directors_correct_order} "
                    f"(-CH3 < -C2H5 < -Cl). The proposed order for these was {llm_op_order}.")

        # Generic failure message
        return f"The proposed order {llm_order} is incorrect. The correct order is {expected_order}."

# Execute the check
result = check_electrophilic_substitution_order()
print(result)