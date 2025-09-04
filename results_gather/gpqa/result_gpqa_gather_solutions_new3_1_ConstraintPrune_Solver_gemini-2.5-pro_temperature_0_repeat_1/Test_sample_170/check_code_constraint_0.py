import re

def check_answer():
    """
    Checks the correctness of the answer for the electrophilic substitution question.

    The function encodes the chemical principles for determining the para-isomer yield:
    1.  Substituents are classified as meta-directing or ortho,para-directing.
    2.  Meta-directors always yield less para-isomer than ortho,para-directors.
    3.  Within meta-directors, yield is inversely related to deactivating strength.
    4.  Within ortho,para-directors, yield is related to steric hindrance and electronic effects favoring the para position.
    """
    
    # The final answer from the LLM to be checked.
    llm_answer_option = "D"
    llm_answer_sequence_str = "4 < 6 < 2 < 1 < 5 < 3"

    # --- Step 1: Define the chemical properties and rules ---
    # Lower 'rank' means lower para-isomer yield within the group.
    substance_properties = {
        # Meta-directors (lowest para-yield)
        # Order is based on the reverse of deactivating strength: -NO2 > -COOH > -COOC2H5
        4: {'name': 'Nitrobenzene (-NO2)', 'director': 'meta', 'rank_within_group': 1},
        6: {'name': 'Benzoic acid (-COOH)', 'director': 'meta', 'rank_within_group': 2},
        2: {'name': 'Ethyl benzoate (-COOC2H5)', 'director': 'meta', 'rank_within_group': 3},

        # Ortho, Para-directors (highest para-yield)
        # Order is based on para-selectivity (sterics and electronics)
        1: {'name': 'Toluene (-CH3)', 'director': 'op', 'rank_within_group': 1},
        5: {'name': 'Ethylbenzene (-C2H5)', 'director': 'op', 'rank_within_group': 2},
        3: {'name': 'Chlorobenzene (-Cl)', 'director': 'op', 'rank_within_group': 3},
    }

    # --- Step 2: Generate the correct sequence based on the rules ---
    meta_directors = []
    op_directors = []

    for sub_id, props in substance_properties.items():
        if props['director'] == 'meta':
            meta_directors.append((sub_id, props))
        else:
            op_directors.append((sub_id, props))

    # Sort each group based on their rank
    sorted_meta = sorted(meta_directors, key=lambda x: x[1]['rank_within_group'])
    sorted_op = sorted(op_directors, key=lambda x: x[1]['rank_within_group'])

    # Extract the IDs to form the correct sequence
    correct_sequence_ids = [item[0] for item in sorted_meta] + [item[0] for item in sorted_op]

    # --- Step 3: Parse the LLM's answer and compare ---
    try:
        # Extract numbers from the string like "4 < 6 < 2 < 1 < 5 < 3"
        llm_sequence_ids = [int(s) for s in re.findall(r'\d+', llm_answer_sequence_str)]
    except (ValueError, IndexError):
        return f"Could not parse the sequence from the provided answer: '{llm_answer_sequence_str}'"

    if len(llm_sequence_ids) != 6:
        return f"The provided answer sequence does not contain 6 substances. Found: {llm_sequence_ids}"

    # --- Step 4: Check for correctness and provide reasons if incorrect ---
    if llm_sequence_ids == correct_sequence_ids:
        return "Correct"
    else:
        # Check Constraint 1: Meta-directors must come before ortho,para-directors
        llm_meta_part = llm_sequence_ids[:3]
        llm_op_part = llm_sequence_ids[3:]
        
        correct_meta_ids = {item[0] for item in sorted_meta}
        correct_op_ids = {item[0] for item in sorted_op}

        if not set(llm_meta_part) == correct_meta_ids:
            incorrect_placements = [substance_properties[i]['name'] for i in llm_meta_part if i not in correct_meta_ids]
            return (f"Incorrect. The fundamental grouping is wrong. "
                    f"The first three substances should be the meta-directors (Nitrobenzene, Benzoic acid, Ethyl benzoate), "
                    f"but the answer places {', '.join(incorrect_placements)} in the first half of the sequence.")

        # Check Constraint 2: Ordering within meta-directors
        correct_meta_order = [item[0] for item in sorted_meta]
        if llm_meta_part != correct_meta_order:
            return (f"Incorrect. The relative ordering of the meta-directing groups is wrong. "
                    f"The order of increasing para-yield is determined by the reverse of deactivating strength (-NO2 < -COOH < -COOC2H5). "
                    f"The correct order for this group is {correct_meta_order}, but the answer provides {llm_meta_part}.")

        # Check Constraint 3: Ordering within ortho,para-directors
        correct_op_order = [item[0] for item in sorted_op]
        if llm_op_part != correct_op_order:
            return (f"Incorrect. The relative ordering of the ortho,para-directing groups is wrong. "
                    f"The order of increasing para-yield is Toluene < Ethylbenzene < Chlorobenzene due to steric and electronic effects. "
                    f"The correct order for this group is {correct_op_order}, but the answer provides {llm_op_part}.")
        
        return "Incorrect. The overall sequence is wrong for reasons not captured by the specific checks."


# Run the check
result = check_answer()
print(result)