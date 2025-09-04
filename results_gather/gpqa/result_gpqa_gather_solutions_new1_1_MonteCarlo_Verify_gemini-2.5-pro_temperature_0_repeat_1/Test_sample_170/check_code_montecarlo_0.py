def check_chemistry_answer():
    """
    This function checks the correctness of the provided answer for a chemistry question.
    The question asks to arrange six substituted benzenes in order of increasing
    weight fraction of the yield of the para-isomer upon electrophilic bromination.

    The provided answer to check is C, which corresponds to the sequence: 4 < 6 < 2 < 1 < 5 < 3.
    """

    # The sequence from the provided answer 'C'
    answer_sequence = [4, 6, 2, 1, 5, 3]

    # --- Constraint 1: Correct Grouping by Directing Effect ---
    # The sequence must be divided into meta-directors (lowest para-yield) followed by
    # ortho/para-directors (highest para-yield).
    meta_directors_set = {2, 4, 6}
    op_directors_set = {1, 3, 5}

    first_three_in_answer = set(answer_sequence[:3])
    last_three_in_answer = set(answer_sequence[3:])

    if not (first_three_in_answer == meta_directors_set and last_three_in_answer == op_directors_set):
        return (f"Constraint Violated: Incorrect grouping. "
                f"The substances with meta-directing groups ({meta_directors_set}) should all come before "
                f"the substances with ortho/para-directing groups ({op_directors_set}). "
                f"The answer places {first_three_in_answer} first, which is incorrect.")

    # --- Constraint 2: Order of Meta-Directors ---
    # For meta-directors, the para-yield is inversely related to the deactivating strength.
    # Order of deactivating strength: -NO2 (4) > -COOH (6) > -COOC2H5 (2).
    # Therefore, the order of increasing para-yield must be: 4 < 6 < 2.
    meta_part_in_answer = answer_sequence[:3]
    correct_meta_order = [4, 6, 2]

    if meta_part_in_answer != correct_meta_order:
        return (f"Constraint Violated: Incorrect order of meta-directors. "
                f"Based on deactivating strength (-NO2 > -COOH > -COOC2H5), the order of increasing "
                f"para-isomer yield should be {correct_meta_order}. "
                f"The answer provides {meta_part_in_answer}.")

    # --- Constraint 3: Order of Ortho/Para-Directors ---
    # For ortho/para-directors, the para-yield is determined by steric and electronic effects.
    # Steric hindrance: -C2H5 is bulkier than -CH3, so para-yield(5) > para-yield(1).
    # Halogen selectivity: -Cl is highly para-selective, more so than alkyl groups.
    # Therefore, the order of increasing para-yield must be: 1 < 5 < 3.
    op_part_in_answer = answer_sequence[3:]
    correct_op_order = [1, 5, 3]

    if op_part_in_answer != correct_op_order:
        return (f"Constraint Violated: Incorrect order of ortho/para-directors. "
                f"Based on steric hindrance and electronic effects, the order of increasing "
                f"para-isomer yield should be {correct_op_order} (Toluene < Ethylbenzene < Chlorobenzene). "
                f"The answer provides {op_part_in_answer}.")

    # If all constraints are satisfied, the answer is correct.
    return "Correct"

# Execute the check
result = check_chemistry_answer()
print(result)