import re

def check_answer():
    """
    Checks the correctness of the answer for the given chemistry question.

    The question asks to arrange substances in order of increasing weight fraction
    of the para-isomer yield during electrophilic bromination.

    The logic is based on established principles of electrophilic aromatic substitution:
    1.  Substituents are classified as meta-directing or ortho,para-directing.
    2.  Meta-directors yield very little para-product and are ranked first.
    3.  Ortho,para-directors yield significant para-product and are ranked last.
    4.  Within meta-directors, the para-yield is inversely proportional to the
        group's deactivating strength (-NO2 > -COOH > -COOC2H5).
    5.  Within ortho,para-directors, the para-yield is influenced by steric
        and electronic effects, generally increasing with steric bulk (-CH3 < -C2H5)
        and being highest for halogens (-Cl).
    """
    
    # Define the substances and their properties relevant to the question
    # Rank is from 1 (lowest para-yield) to 6 (highest para-yield)
    substances = {
        1: {'name': 'Toluene', 'substituent': '-CH3', 'type': 'o,p', 'rank': 4},
        2: {'name': 'Ethyl benzoate', 'substituent': '-COOC2H5', 'type': 'meta', 'rank': 3},
        3: {'name': 'Chlorobenzene', 'substituent': '-Cl', 'type': 'o,p', 'rank': 6},
        4: {'name': 'Nitrobenzene', 'substituent': '-NO2', 'type': 'meta', 'rank': 1},
        5: {'name': 'Ethylbenzene', 'substituent': '-C2H5', 'type': 'o,p', 'rank': 5},
        6: {'name': 'Benzoic acid', 'substituent': '-COOH', 'type': 'meta', 'rank': 2},
    }

    # The correct order based on chemical principles
    correct_order = sorted(substances.keys(), key=lambda k: substances[k]['rank'])
    
    # The options provided in the question
    options = {
        'A': [3, 5, 1, 6, 2, 4],
        'B': [4, 6, 2, 1, 5, 3],
        'C': [4, 2, 6, 3, 1, 5],
        'D': [6, 2, 4, 5, 1, 3]
    }

    # The candidate answer to check
    candidate_answer_letter = 'B'
    candidate_order = options.get(candidate_answer_letter)

    if not candidate_order:
        return f"Invalid option '{candidate_answer_letter}'. Please choose from A, B, C, D."

    # Check if the candidate order matches the correct order
    if candidate_order == correct_order:
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        
        # Check the fundamental grouping
        meta_directors = {k for k, v in substances.items() if v['type'] == 'meta'}
        op_directors = {k for k, v in substances.items() if v['type'] == 'o,p'}
        
        candidate_meta_part = [x for x in candidate_order if x in meta_directors]
        candidate_op_part = [x for x in candidate_order if x in op_directors]
        
        if candidate_order[:3] != candidate_meta_part or candidate_order[3:] != candidate_op_part:
            return (f"Incorrect. The fundamental grouping is wrong. "
                    f"Meta-directing substances (which give low para-yield) should come before "
                    f"ortho,para-directing substances (which give high para-yield). "
                    f"The correct grouping is {list(meta_directors)} followed by {list(op_directors)}.")

        # Check the order within meta-directors
        correct_meta_order = sorted(meta_directors, key=lambda k: substances[k]['rank'])
        if candidate_meta_part != correct_meta_order:
            return (f"Incorrect. The relative ordering of the meta-directing substances is wrong. "
                    f"The para-yield is inversely related to the deactivating strength of the group. "
                    f"The order of deactivating strength is -NO2 (4) > -COOH (6) > -COOC2H5 (2). "
                    f"Therefore, the correct order of increasing para-yield is {correct_meta_order}, "
                    f"but the answer provided {candidate_meta_part}.")

        # Check the order within ortho,para-directors
        correct_op_order = sorted(op_directors, key=lambda k: substances[k]['rank'])
        if candidate_op_part != correct_op_order:
            return (f"Incorrect. The relative ordering of the ortho,para-directing substances is wrong. "
                    f"Para-selectivity increases with the substituent's steric bulk and is very high for halogens. "
                    f"The correct order of increasing para-yield is Toluene (1) < Ethylbenzene (5) < Chlorobenzene (3). "
                    f"Therefore, the correct order is {correct_op_order}, "
                    f"but the answer provided {candidate_op_part}.")
        
        return f"Incorrect. The final proposed order {candidate_order} does not match the correct order {correct_order}."

# Run the check
result = check_answer()
print(result)