def check_chemistry_answer():
    """
    Checks the correctness of a proposed order for para-isomer yield in electrophilic bromination.
    The function encodes chemical principles and compares the given answer against them.
    """
    
    # --- 1. Define Chemical Data and Principles ---

    # Information about each substituent on the benzene ring.
    # 'type' classifies the substituent's directing effect.
    substance_properties = {
        1: {'name': 'Toluene (C6H5-CH3)', 'type': 'ortho_para'},
        2: {'name': 'Ethyl benzoate (C6H5-COOC2H5)', 'type': 'meta'},
        3: {'name': 'Chlorobenzene (C6H5-Cl)', 'type': 'ortho_para'},
        4: {'name': 'Nitrobenzene (C6H5-NO2)', 'type': 'meta'},
        5: {'name': 'Ethylbenzene (C6H5-C2H5)', 'type': 'ortho_para'},
        6: {'name': 'Benzoic acid (C6H5-COOH)', 'type': 'meta'},
    }

    # The answer to be checked, corresponding to option A: 4 < 6 < 2 < 1 < 5 < 3
    proposed_order = (4, 6, 2, 1, 5, 3)

    # --- 2. Define Correct Ordering Based on Chemical Rules ---

    # Rule A: Meta-directing groups yield very little para-isomer compared to ortho/para-directing groups.
    # Therefore, all meta-directors must come before all ortho/para-directors in the final order.
    meta_directors = {k for k, v in substance_properties.items() if v['type'] == 'meta'}
    op_directors = {k for k, v in substance_properties.items() if v['type'] == 'ortho_para'}

    # Rule B: Ordering within the meta-directing group.
    # The para-yield is inversely related to the deactivating strength of the group.
    # Deactivating strength: -NO2 > -COOH > -COOC2H5
    # Expected increasing para-yield order: Nitrobenzene (4) < Benzoic acid (6) < Ethyl benzoate (2)
    expected_meta_order = (4, 6, 2)

    # Rule C: Ordering within the ortho/para-directing group.
    # The order depends on the para/ortho ratio, influenced by sterics and electronics.
    # - Steric Hindrance: The bulkier ethyl group (-C2H5) hinders ortho-attack more than the methyl group (-CH3),
    #   increasing the relative yield of the para-isomer. So, para(Toluene) < para(Ethylbenzene). Order: 1 < 5.
    # - Electronic Effects: The -Cl group, while deactivating overall, strongly favors the para position for
    #   electrophilic attack, resulting in a very high para-isomer yield (~89%), typically higher than for alkylbenzenes.
    # Expected increasing para-yield order: Toluene (1) < Ethylbenzene (5) < Chlorobenzene (3)
    expected_op_order = (1, 5, 3)

    # --- 3. Perform the Checks ---
    
    errors = []

    # Check 1: Grouping of meta vs. ortho/para directors
    num_meta = len(meta_directors)
    found_meta_group = set(proposed_order[:num_meta])
    found_op_group = set(proposed_order[num_meta:])

    if found_meta_group != meta_directors:
        errors.append(f"Constraint Failure: Grouping. The first {num_meta} substances should be the meta-directors {sorted(list(meta_directors))}, but the answer proposes {sorted(list(found_meta_group))}.")
    
    if found_op_group != op_directors:
        errors.append(f"Constraint Failure: Grouping. The last {len(op_directors)} substances should be the ortho/para-directors {sorted(list(op_directors))}, but the answer proposes {sorted(list(found_op_group))}.")

    # Check 2: Internal order of meta-directors
    found_meta_order = proposed_order[:num_meta]
    if found_meta_order != expected_meta_order:
        errors.append(f"Constraint Failure: Meta-director order. Expected {expected_meta_order} based on deactivating strength, but the answer proposes {found_meta_order}.")

    # Check 3: Internal order of ortho/para-directors
    found_op_order = proposed_order[num_meta:]
    if found_op_order != expected_op_order:
        errors.append(f"Constraint Failure: Ortho/para-director order. Expected {expected_op_order} based on steric/electronic effects, but the answer proposes {found_op_order}.")

    # --- 4. Return Final Verdict ---
    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check and print the result
result = check_chemistry_answer()
print(result)