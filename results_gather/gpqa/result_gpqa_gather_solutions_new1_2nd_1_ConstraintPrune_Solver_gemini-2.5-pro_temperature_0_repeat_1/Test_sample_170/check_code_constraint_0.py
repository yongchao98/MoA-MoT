import re

def check_answer_correctness(llm_answer_text):
    """
    Checks the correctness of the LLM's answer for the given chemistry question.

    The function verifies the answer based on established principles of electrophilic
    aromatic substitution:
    1.  Separation of substituents into meta-directors and ortho,para-directors.
    2.  Ordering within meta-directors based on deactivating strength.
    3.  Ordering within ortho,para-directors based on steric and electronic effects.
    """

    # Define the chemical properties of each substance
    substances = {
        # id: {'name': str, 'directing_effect': 'meta' or 'op'}
        1: {'name': 'Toluene (-CH3)', 'directing_effect': 'op'},
        2: {'name': 'Ethyl benzoate (-COOC2H5)', 'directing_effect': 'meta'},
        3: {'name': 'Chlorobenzene (-Cl)', 'directing_effect': 'op'},
        4: {'name': 'Nitrobenzene (-NO2)', 'directing_effect': 'meta'},
        5: {'name': 'Ethylbenzene (-C2H5)', 'directing_effect': 'op'},
        6: {'name': 'Benzoic acid (-COOH)', 'directing_effect': 'meta'},
    }

    # Define the correct ordering based on chemical principles
    # Meta-directors (low para-yield):
    # Deactivating strength: -NO2 > -COOH > -COOC2H5
    # Inverse order for para-yield: 4 < 6 < 2
    correct_meta_order = [4, 6, 2]

    # Ortho,para-directors (high para-yield):
    # Steric hindrance: -C2H5 > -CH3 => para-yield(5) > para-yield(1)
    # Halogen effect: -Cl is highly para-selective => para-yield(3) is highest
    # Order for para-yield: 1 < 5 < 3
    correct_op_order = [1, 5, 3]

    # The full correct sequence
    correct_full_order = correct_meta_order + correct_op_order

    # Map options to their sequences
    options = {
        "A": [3, 5, 1, 6, 2, 4],
        "B": [4, 6, 2, 1, 5, 3],
        "C": [6, 2, 4, 5, 1, 3],
        "D": [4, 2, 6, 3, 1, 5]
    }

    # --- Verification Logic ---
    # 1. Extract the chosen option key (e.g., 'B') from the LLM's response
    match = re.search(r'<<<([A-D])>>>', llm_answer_text)
    if not match:
        return "Incorrect. The answer is not in the required format '<<<X>>>' where X is A, B, C, or D."
    
    answer_key = match.group(1)
    proposed_order = options[answer_key]

    # 2. Check if the proposed order matches the correct one
    if proposed_order == correct_full_order:
        return "Correct"

    # 3. If incorrect, provide a specific reason
    # Constraint 1: Meta-directors must all come before ortho,para-directors
    meta_ids = {k for k, v in substances.items() if v['directing_effect'] == 'meta'}
    op_ids = {k for k, v in substances.items() if v['directing_effect'] == 'op'}
    
    last_meta_index = -1
    first_op_index = len(proposed_order)
    for i, sub_id in enumerate(proposed_order):
        if sub_id in meta_ids:
            last_meta_index = i
        if sub_id in op_ids and i < first_op_index:
            first_op_index = i
            
    if last_meta_index > first_op_index:
        return (f"Incorrect. The answer '{answer_key}' proposes the order {proposed_order}. "
                f"This violates the primary rule that all meta-directing compounds ({meta_ids}) must have lower "
                f"para-yields than all ortho,para-directing compounds ({op_ids}).")

    # Constraint 2: Check the order within the meta-directing group
    proposed_meta_order = [sub_id for sub_id in proposed_order if sub_id in meta_ids]
    if proposed_meta_order != correct_meta_order:
        return (f"Incorrect. The answer '{answer_key}' proposes the order {proposed_order}. "
                f"The relative ordering of the meta-directing compounds is wrong. "
                f"Para-yield is inversely related to deactivating strength (-NO2 > -COOH > -COOC2H5), "
                f"so the correct sub-order is {correct_meta_order}. The proposed sub-order is {proposed_meta_order}.")

    # Constraint 3: Check the order within the ortho,para-directing group
    proposed_op_order = [sub_id for sub_id in proposed_order if sub_id in op_ids]
    if proposed_op_order != correct_op_order:
        return (f"Incorrect. The answer '{answer_key}' proposes the order {proposed_order}. "
                f"The relative ordering of the ortho,para-directing compounds is wrong. "
                f"Based on steric hindrance and electronic effects, the correct sub-order is {correct_op_order} "
                f"(Toluene < Ethylbenzene < Chlorobenzene). The proposed sub-order is {proposed_op_order}.")

    return "Incorrect. The proposed order is wrong for reasons not explicitly checked."

# The final answer provided by the LLM to be checked
final_answer_from_prompt = "<<<B>>>"

# Run the check
result = check_answer_correctness(final_answer_from_prompt)
print(result)