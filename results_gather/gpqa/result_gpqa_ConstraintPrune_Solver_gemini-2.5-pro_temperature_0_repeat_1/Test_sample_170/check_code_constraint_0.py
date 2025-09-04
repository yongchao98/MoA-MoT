def check_chemistry_answer():
    """
    Checks the correctness of the proposed order for electrophilic bromination para-isomer yield.

    The function models the chemical principles as a set of constraints:
    1.  Substances with meta-directing groups have lower para-yield than those with ortho,para-directing groups.
    2.  Within meta-directors, the para-yield is inversely related to the group's deactivating strength.
        (Deactivating strength: NO2 > COOH > COOC2H5 => Para-yield: NO2 < COOH < COOC2H5)
    3.  Within ortho,para-directors, the para-yield is determined by steric and electronic effects.
        (Para-selectivity: Cl > C2H5 > CH3)
    """
    # Define substance properties. Ranks indicate the expected position in the final sorted list within their group.
    # Lower rank = lower para-yield.
    substances = {
        1: {'name': 'Toluene', 'substituent': '-CH3', 'director': 'op', 'rank': 1},
        2: {'name': 'Ethyl benzoate', 'substituent': '-COOC2H5', 'director': 'meta', 'rank': 3},
        3: {'name': 'Chlorobenzene', 'substituent': '-Cl', 'director': 'op', 'rank': 3},
        4: {'name': 'Nitrobenzene', 'substituent': '-NO2', 'director': 'meta', 'rank': 1},
        5: {'name': 'Ethylbenzene', 'substituent': '-C2H5', 'director': 'op', 'rank': 2},
        6: {'name': 'Benzoic acid', 'substituent': '-COOH', 'director': 'meta', 'rank': 2},
    }

    # The proposed answer from the LLM corresponds to option B.
    proposed_answer = [4, 6, 2, 1, 5, 3]

    # Constraint 1: All meta-directors must appear before all ortho,para-directors.
    meta_group_indices = []
    op_group_indices = []
    for i, sub_id in enumerate(proposed_answer):
        if substances[sub_id]['director'] == 'meta':
            meta_group_indices.append(i)
        else:
            op_group_indices.append(i)

    if not meta_group_indices or not op_group_indices:
         return "Constraint check failed: The answer must contain both meta and ortho,para directing groups."

    if max(meta_group_indices) > min(op_group_indices):
        return (f"Constraint check failed: Grouping is incorrect. "
                f"The last meta-director (at index {max(meta_group_indices)}) appears after "
                f"the first ortho,para-director (at index {min(op_group_indices)}).")

    # Constraint 2: Check the internal order of the meta-directing group.
    meta_sub_ids = [sub_id for sub_id in proposed_answer if substances[sub_id]['director'] == 'meta']
    expected_meta_order = sorted(meta_sub_ids, key=lambda sub_id: substances[sub_id]['rank'])
    if meta_sub_ids != expected_meta_order:
        actual_str = ' < '.join(map(str, meta_sub_ids))
        expected_str = ' < '.join(map(str, expected_meta_order))
        return (f"Constraint check failed: The order within the meta-directing substances is incorrect. "
                f"Expected order based on increasing para-yield is {expected_str}, but the answer gives {actual_str}.")

    # Constraint 3: Check the internal order of the ortho,para-directing group.
    op_sub_ids = [sub_id for sub_id in proposed_answer if substances[sub_id]['director'] == 'op']
    expected_op_order = sorted(op_sub_ids, key=lambda sub_id: substances[sub_id]['rank'])
    if op_sub_ids != expected_op_order:
        actual_str = ' < '.join(map(str, op_sub_ids))
        expected_str = ' < '.join(map(str, expected_op_order))
        return (f"Constraint check failed: The order within the ortho,para-directing substances is incorrect. "
                f"Expected order based on increasing para-yield is {expected_str}, but the answer gives {actual_str}.")

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check
result = check_chemistry_answer()
# This will print "Correct" if the logic holds, or an error message otherwise.
# print(result) # You can uncomment this to see the output when running the script.