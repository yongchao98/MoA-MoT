import collections

def check_answer():
    """
    Checks the correctness of the final answer for the electrophilic substitution question.

    The function encodes the chemical principles for determining the yield of the para-isomer
    and compares the resulting correct order with the provided answer.
    """

    # Define the substances and their properties based on chemical principles.
    # 'director_type': 'meta' or 'op' (ortho,para). Meta-directors have lower para-yield.
    # 'rank_in_type': A rank within each director type. Lower rank means lower para-yield.
    substances = [
        # --- Meta-Directors ---
        # The order of increasing para-yield is inversely related to deactivating strength.
        # Deactivating strength: -NO2 > -COOH > -COOC2H5
        # Therefore, para-yield order: NO2 < COOH < COOC2H5
        {'id': 4, 'name': 'Nitrobenzene', 'substituent': '-NO2', 'director_type': 'meta', 'rank_in_type': 1},
        {'id': 6, 'name': 'Benzoic acid', 'substituent': '-COOH', 'director_type': 'meta', 'rank_in_type': 2},
        {'id': 2, 'name': 'Ethyl benzoate', 'substituent': '-COOC2H5', 'director_type': 'meta', 'rank_in_type': 3},

        # --- Ortho, Para-Directors ---
        # The order of increasing para-yield is based on steric hindrance and electronic effects.
        # Para-selectivity: -Cl > -C2H5 > -CH3
        # Therefore, para-yield order: CH3 < C2H5 < Cl
        {'id': 1, 'name': 'Toluene', 'substituent': '-CH3', 'director_type': 'op', 'rank_in_type': 1},
        {'id': 5, 'name': 'Ethylbenzene', 'substituent': '-C2H5', 'director_type': 'op', 'rank_in_type': 2},
        {'id': 3, 'name': 'Chlorobenzene', 'substituent': '-Cl', 'director_type': 'op', 'rank_in_type': 3},
    ]

    # Sort the substances based on the encoded chemical principles to get the correct order.
    # First, sort by director type ('meta' comes before 'op').
    # Then, sort by the rank within that type.
    correctly_sorted = sorted(substances, key=lambda x: (x['director_type'], x['rank_in_type']))
    correct_sequence = [s['id'] for s in correctly_sorted]

    # The final answer provided was 'C', which corresponds to the sequence 4 < 6 < 2 < 1 < 5 < 3.
    # Let's define this sequence for checking.
    answer_sequence = [4, 6, 2, 1, 5, 3]

    # --- Verification ---
    errors = []

    # Constraint 1: Meta-directors must come before ortho,para-directors.
    # The first 3 in the sequence should be meta-directors, the last 3 should be o,p-directors.
    meta_ids = {s['id'] for s in substances if s['director_type'] == 'meta'}
    op_ids = {s['id'] for s in substances if s['director_type'] == 'op'}
    
    if set(answer_sequence[:3]) != meta_ids:
        errors.append(f"Constraint Violation: The first three substances should be the meta-directors {sorted(list(meta_ids))}, but the answer provides {answer_sequence[:3]}.")
    
    if set(answer_sequence[3:]) != op_ids:
        errors.append(f"Constraint Violation: The last three substances should be the ortho,para-directors {sorted(list(op_ids))}, but the answer provides {answer_sequence[3:]}.")

    # Constraint 2: The internal order of meta-directors must be correct.
    # Expected order: 4 < 6 < 2
    correct_meta_order = [s['id'] for s in correctly_sorted if s['director_type'] == 'meta']
    if answer_sequence[:3] != correct_meta_order:
        errors.append(f"Constraint Violation: The order of meta-directors is incorrect. Based on deactivating strength, the order of increasing para-yield should be {correct_meta_order}, but the answer provides {answer_sequence[:3]}.")

    # Constraint 3: The internal order of ortho,para-directors must be correct.
    # Expected order: 1 < 5 < 3
    correct_op_order = [s['id'] for s in correctly_sorted if s['director_type'] == 'op']
    if answer_sequence[3:] != correct_op_order:
        errors.append(f"Constraint Violation: The order of ortho,para-directors is incorrect. Based on steric/electronic effects, the order of increasing para-yield should be {correct_op_order}, but the answer provides {answer_sequence[3:]}.")

    # Final check on the whole sequence
    if correct_sequence != answer_sequence and not errors:
        # This is a fallback, the specific checks above should catch any errors.
        errors.append(f"The final sequence is incorrect. Expected {correct_sequence} but got {answer_sequence}.")

    if not errors:
        return "Correct"
    else:
        return "\n".join(errors)

# Run the check and print the result.
result = check_answer()
print(result)