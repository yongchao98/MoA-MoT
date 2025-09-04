def check_electrophilic_substitution_order():
    """
    Checks the correctness of the proposed order for para-isomer yield in electrophilic bromination.

    The function verifies the answer against three core chemical principles:
    1. Meta-directing compounds yield less para-isomer than ortho,para-directing compounds.
    2. The order among meta-directors is determined by their deactivating strength.
    3. The order among ortho,para-directors is determined by steric and electronic effects influencing para-selectivity.
    """
    # Define the properties of each substance.
    # 'type': 'meta' or 'op' (ortho,para)
    # 'rank': A numerical rank within each type, where a lower number means a lower para-yield.
    compounds = {
        1: {'name': 'Toluene', 'substituent': '-CH3', 'type': 'op', 'rank': 1},
        2: {'name': 'Ethyl benzoate', 'substituent': '-COOC2H5', 'type': 'meta', 'rank': 3},
        3: {'name': 'Chlorobenzene', 'substituent': '-Cl', 'type': 'op', 'rank': 3},
        4: {'name': 'Nitrobenzene', 'substituent': '-NO2', 'type': 'meta', 'rank': 1},
        5: {'name': 'Ethylbenzene', 'substituent': '-C2H5', 'type': 'op', 'rank': 2},
        6: {'name': 'Benzoic acid', 'substituent': '-COOH', 'type': 'meta', 'rank': 2},
    }

    # The proposed answer from the LLM (Option D)
    proposed_order = [4, 6, 2, 1, 5, 3]

    # --- Verification Step 1: Grouping by Directing Effect ---
    meta_directors_in_order = [c for c in proposed_order if compounds[c]['type'] == 'meta']
    op_directors_in_order = [c for c in proposed_order if compounds[c]['type'] == 'op']

    if proposed_order != meta_directors_in_order + op_directors_in_order:
        return ("Incorrect: The order does not correctly group the compounds. "
                "All meta-directing compounds (4, 6, 2) must be placed before all "
                "ortho,para-directing compounds (1, 5, 3) because they produce a "
                "significantly lower yield of the para-isomer.")

    # --- Verification Step 2: Ordering of Meta-Directors ---
    # The order should follow the ranks: 1 < 2 < 3
    # This corresponds to deactivating strength: -NO2 > -COOH > -COOC2H5
    actual_meta_ranks = [compounds[c]['rank'] for c in meta_directors_in_order]
    if actual_meta_ranks != sorted(actual_meta_ranks):
        correct_meta_order = sorted(meta_directors_in_order, key=lambda c: compounds[c]['rank'])
        return (f"Incorrect: The relative order of the meta-directing compounds is wrong. "
                f"Based on deactivating strength, the correct order of increasing para-yield is {correct_meta_order}, "
                f"but the answer provides {meta_directors_in_order}.")

    # --- Verification Step 3: Ordering of Ortho,Para-Directors ---
    # The order should follow the ranks: 1 < 2 < 3
    # This corresponds to para-selectivity: Toluene < Ethylbenzene < Chlorobenzene
    actual_op_ranks = [compounds[c]['rank'] for c in op_directors_in_order]
    if actual_op_ranks != sorted(actual_op_ranks):
        correct_op_order = sorted(op_directors_in_order, key=lambda c: compounds[c]['rank'])
        return (f"Incorrect: The relative order of the ortho,para-directing compounds is wrong. "
                f"Based on steric and electronic effects, the correct order of increasing para-yield is {correct_op_order}, "
                f"but the answer provides {op_directors_in_order}.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_electrophilic_substitution_order()
print(result)