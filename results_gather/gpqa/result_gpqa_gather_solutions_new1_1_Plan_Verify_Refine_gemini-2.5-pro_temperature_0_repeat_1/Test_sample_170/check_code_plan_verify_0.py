import collections

def check_answer():
    """
    Checks the correctness of the provided answer for the chemistry question.

    The question asks to arrange substances in order of increasing weight fraction
    of the yield of the para-isomer in an electrophilic bromination reaction.

    The substances are:
    1) C6H5-CH3 (Toluene)
    2) C6H5-COOC2H5 (Ethyl benzoate)
    3) C6H5-Cl (Chlorobenzene)
    4) C6H5-NO2 (Nitrobenzene)
    5) C6H5-C2H5 (Ethylbenzene)
    6) C6H5-COOH (Benzoic acid)

    The provided answer is D, which corresponds to the sequence 4 < 6 < 2 < 1 < 5 < 3.
    """

    # The provided answer from the LLM to be checked
    llm_answer_option = 'D'

    # The options given in the question
    options = {
        'A': [3, 5, 1, 6, 2, 4],
        'B': [4, 2, 6, 3, 1, 5],
        'C': [6, 2, 4, 5, 1, 3],
        'D': [4, 6, 2, 1, 5, 3]
    }

    if llm_answer_option not in options:
        return f"Invalid option '{llm_answer_option}' provided. The option must be one of {list(options.keys())}."

    candidate_order = options[llm_answer_option]

    # --- Verification based on chemical principles ---

    # 1. Define the properties of each substituent
    # 'group' classifies them as ortho/para or meta directors.
    substances = {
        1: {'name': 'Toluene', 'group': 'op_director'},
        2: {'name': 'Ethyl benzoate', 'group': 'meta_director'},
        3: {'name': 'Chlorobenzene', 'group': 'op_director'},
        4: {'name': 'Nitrobenzene', 'group': 'meta_director'},
        5: {'name': 'Ethylbenzene', 'group': 'op_director'},
        6: {'name': 'Benzoic acid', 'group': 'meta_director'}
    }
    meta_directors = {k for k, v in substances.items() if v['group'] == 'meta_director'}
    op_directors = {k for k, v in substances.items() if v['group'] == 'op_director'}

    # Constraint 1: All meta-directors must have a lower para-yield than all ortho/para-directors.
    # This means the first three substances in the order must be the meta-directors.
    first_three = set(candidate_order[:3])
    last_three = set(candidate_order[3:])
    
    if not (first_three == meta_directors and last_three == op_directors):
        return (f"Incorrect. Constraint violated: The order must group all meta-directing substances before all ortho/para-directing substances. "
                f"Meta-directors are {{2, 4, 6}} and ortho/para-directors are {{1, 3, 5}}. "
                f"The given order {candidate_order} does not follow this grouping.")

    # Constraint 2: Check the relative order of meta-directors.
    # The para-yield is inversely related to the deactivating strength of the group.
    # Deactivating strength: -NO2 > -COOH > -COOC2H5
    # Therefore, increasing para-yield order: 4 < 6 < 2
    meta_sub_order = [s for s in candidate_order if s in meta_directors]
    correct_meta_order = [4, 6, 2]
    if meta_sub_order != correct_meta_order:
        return (f"Incorrect. Constraint violated: The relative order of meta-directing substances is wrong. "
                f"The expected order based on increasing para-yield is {correct_meta_order} (Nitrobenzene < Benzoic acid < Ethyl benzoate). "
                f"The given order has them as {meta_sub_order}.")

    # Constraint 3: Check the relative order of ortho/para-directors.
    # The para-yield is influenced by steric hindrance and electronic effects.
    # Sterics: -C2H5 is bulkier than -CH3, so para-yield(Toluene) < para-yield(Ethylbenzene). (1 < 5)
    # Electronics/Data: -Cl gives a very high para-yield, higher than alkylbenzenes.
    # Therefore, increasing para-yield order: 1 < 5 < 3
    op_sub_order = [s for s in candidate_order if s in op_directors]
    correct_op_order = [1, 5, 3]
    if op_sub_order != correct_op_order:
        return (f"Incorrect. Constraint violated: The relative order of ortho/para-directing substances is wrong. "
                f"The expected order based on increasing para-yield is {correct_op_order} (Toluene < Ethylbenzene < Chlorobenzene). "
                f"The given order has them as {op_sub_order}.")

    # If all constraints are passed, the answer is correct.
    return "Correct"

# Run the check
result = check_answer()
print(result)