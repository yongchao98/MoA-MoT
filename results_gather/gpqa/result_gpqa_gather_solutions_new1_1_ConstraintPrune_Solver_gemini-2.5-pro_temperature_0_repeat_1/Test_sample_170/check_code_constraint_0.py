def check_chemistry_order():
    """
    Checks the correctness of the proposed order for para-isomer yield in electrophilic bromination.

    The function encodes chemical principles as a series of constraints:
    1.  Meta-directing groups yield less para-isomer than ortho/para-directing groups.
    2.  Within meta-directors, the yield is inversely proportional to the deactivating strength.
    3.  Within ortho/para-directors, the yield is influenced by steric hindrance and halogen selectivity.
    """
    # Define the substances and their relevant chemical properties.
    # 'director': 'meta' or 'op' (ortho/para)
    # 'meta_rank': For meta-directors, rank of deactivating strength (1=strongest, thus lowest para-yield).
    # 'op_rank': For o,p-directors, rank of para-yield (1=lowest, 3=highest).
    substances = {
        1: {'name': 'Toluene', 'substituent': '-CH3', 'director': 'op', 'op_rank': 1},
        2: {'name': 'Ethyl benzoate', 'substituent': '-COOC2H5', 'director': 'meta', 'meta_rank': 3},
        3: {'name': 'Chlorobenzene', 'substituent': '-Cl', 'director': 'op', 'op_rank': 3},
        4: {'name': 'Nitrobenzene', 'substituent': '-NO2', 'director': 'meta', 'meta_rank': 1},
        5: {'name': 'Ethylbenzene', 'substituent': '-C2H5', 'director': 'op', 'op_rank': 2},
        6: {'name': 'Benzoic acid', 'substituent': '-COOH', 'director': 'meta', 'meta_rank': 2}
    }

    # The final answer from the LLM is A, which corresponds to the sequence 4 < 6 < 2 < 1 < 5 < 3.
    proposed_sequence = [4, 6, 2, 1, 5, 3]

    # --- Constraint 1: Grouping by Directing Effect ---
    # All meta-directing compounds must appear before all ortho/para-directing compounds.
    meta_indices = [i for i, sub in enumerate(proposed_sequence) if substances[sub]['director'] == 'meta']
    op_indices = [i for i, sub in enumerate(proposed_sequence) if substances[sub]['director'] == 'op']

    if not meta_indices or not op_indices:
        return "Incorrect: The sequence must contain both meta and ortho/para directors."

    if max(meta_indices) > min(op_indices):
        return (f"Incorrect: Constraint Violated (Grouping).\n"
                f"A meta-directing substance (at position {max(meta_indices)+1}) appears after an ortho/para-directing substance (at position {min(op_indices)+1}).\n"
                f"Meta-directors produce significantly less para-isomer and must come first in the sequence.")

    # --- Constraint 2: Order of Meta-Directors ---
    # The order is based on deactivating strength: NO2 > COOH > COOC2H5.
    # This means the para-yield order should be: 4 < 6 < 2.
    meta_subsequence = [sub for sub in proposed_sequence if substances[sub]['director'] == 'meta']
    expected_meta_order = sorted(meta_subsequence, key=lambda s: substances[s]['meta_rank'])

    if meta_subsequence != expected_meta_order:
        return (f"Incorrect: Constraint Violated (Meta-Director Order).\n"
                f"The order of meta-directing substances is incorrect. Based on deactivating strength, the expected order of increasing para-yield is {expected_meta_order}.\n"
                f"The proposed answer has the order {meta_subsequence}.")

    # --- Constraint 3: Order of Ortho/Para-Directors ---
    # The order is based on sterics and halogen selectivity: Toluene < Ethylbenzene < Chlorobenzene.
    # This means the para-yield order should be: 1 < 5 < 3.
    op_subsequence = [sub for sub in proposed_sequence if substances[sub]['director'] == 'op']
    expected_op_order = sorted(op_subsequence, key=lambda s: substances[s]['op_rank'])

    if op_subsequence != expected_op_order:
        return (f"Incorrect: Constraint Violated (Ortho/Para-Director Order).\n"
                f"The order of ortho/para-directing substances is incorrect. Based on steric hindrance and halogen selectivity, the expected order of increasing para-yield is {expected_op_order}.\n"
                f"The proposed answer has the order {op_subsequence}.")

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_order()
print(result)