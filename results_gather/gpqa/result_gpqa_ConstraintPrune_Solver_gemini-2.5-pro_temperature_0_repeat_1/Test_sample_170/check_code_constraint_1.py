def check_electrophilic_substitution_order():
    """
    Checks the correctness of the proposed order for para-isomer yield in bromination.

    The function uses established chemical principles to determine the correct order:
    1.  Groups are classified as meta-directing or ortho,para-directing. Meta-directors
        always yield less para-isomer than ortho,para-directors.
    2.  Within meta-directors, a stronger deactivating group leads to a lower para-yield.
        The deactivating strength is: -NO2 > -COOH > -COOC2H5.
        Therefore, the para-yield order is: -NO2 < -COOH < -COOC2H5.
    3.  Within ortho,para-directors, the para-yield is determined by steric hindrance and
        electronic effects.
        - The ethyl group (-C2H5) is bulkier than the methyl group (-CH3), leading to
          higher para-selectivity for ethylbenzene over toluene.
        - Halogens like chlorine (-Cl) are known to give very high para-selectivity.
        Therefore, the para-yield order is: -CH3 < -C2H5 < -Cl.
    """
    # Define the substances and their properties based on chemical principles.
    # 'type' is the directing effect.
    # 'rank' is a numerical value representing the expected para-yield, where a lower number means less yield.
    substance_data = {
        1: {'name': 'Toluene', 'group': '-CH3', 'type': 'op', 'rank': 4},
        2: {'name': 'Ethyl benzoate', 'group': '-COOC2H5', 'type': 'meta', 'rank': 3},
        3: {'name': 'Chlorobenzene', 'group': '-Cl', 'type': 'op', 'rank': 6},
        4: {'name': 'Nitrobenzene', 'group': '-NO2', 'type': 'meta', 'rank': 1},
        5: {'name': 'Ethylbenzene', 'group': '-C2H5', 'type': 'op', 'rank': 5},
        6: {'name': 'Benzoic acid', 'group': '-COOH', 'type': 'meta', 'rank': 2},
    }

    # The answer to be checked (Option B)
    llm_order = [4, 6, 2, 1, 5, 3]

    # --- Verification ---

    # 1. Generate the correct order based on the defined ranks
    correct_order = sorted(substance_data, key=lambda k: substance_data[k]['rank'])

    if llm_order != correct_order:
        return f"The final order is incorrect. The correct order should be {correct_order}, but the answer provided is {llm_order}."

    # 2. Verify the grouping (meta-directors first, then ortho,para-directors)
    meta_directors = {k for k, v in substance_data.items() if v['type'] == 'meta'}
    op_directors = {k for k, v in substance_data.items() if v['type'] == 'op'}

    llm_meta_part = llm_order[:len(meta_directors)]
    llm_op_part = llm_order[len(meta_directors):]

    if not (set(llm_meta_part) == meta_directors and set(llm_op_part) == op_directors):
        return "Constraint Violated: The answer does not correctly group the substances. All meta-directing substances (4, 6, 2) must come before all ortho,para-directing substances (1, 5, 3)."

    # 3. Verify the order within the meta-directing group
    correct_meta_order = sorted(meta_directors, key=lambda k: substance_data[k]['rank'])
    if llm_meta_part != correct_meta_order:
        return f"Constraint Violated: The relative order of meta-directors is incorrect. The para-yield increases as deactivating strength decreases (-NO2 < -COOH < -COOC2H5). The correct order is {correct_meta_order}, but the answer provided {llm_meta_part}."

    # 4. Verify the order within the ortho,para-directing group
    correct_op_order = sorted(op_directors, key=lambda k: substance_data[k]['rank'])
    if llm_op_part != correct_op_order:
        return f"Constraint Violated: The relative order of ortho,para-directors is incorrect. Based on steric hindrance and selectivity, the correct order is {correct_op_order} (Toluene < Ethylbenzene < Chlorobenzene), but the answer provided {llm_op_part}."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_electrophilic_substitution_order()
print(result)