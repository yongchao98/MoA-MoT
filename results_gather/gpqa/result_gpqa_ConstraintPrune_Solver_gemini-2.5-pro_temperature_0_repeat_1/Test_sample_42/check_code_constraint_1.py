import re

def check_correctness():
    """
    Checks the correctness of the provided answer for the IUPAC naming question.
    """
    # The proposed correct answer from the LLM
    llm_answer_key = "B"

    # All candidate answers
    candidates = {
        "A": "2-(dimethylamino)-3-formyl-4-methoxy-5-cyano-6-hydroxybenzoic acid",
        "B": "5-cyano-2-(dimethylamino)-3-formyl-6-hydroxy-4-methoxybenzoic acid",
        "C": "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid",
        "D": "2-hydroxy-3-cyano-4-methoxy-5-formyl-6-(dimethylamino)benzoic acid"
    }

    # --- Helper functions ---
    def parse_name(name_string):
        """Parses an IUPAC name into a dictionary of {position: substituent}."""
        if "benzoic acid" not in name_string:
            return None, "Name does not contain 'benzoic acid' parent."
        
        substituents = {1: "COOH"}
        substituent_part = name_string.split('benzoic acid')[0]
        
        # Regex to find number-substituent pairs. Handles groups like (dimethylamino).
        # It also captures the full text for alphabetical order checking.
        pattern = r'(\d+)-(\(?[a-zA-Z]+\)?)'
        matches = re.findall(pattern, substituent_part)
        
        parsed_substituents = []
        for pos, group_name in matches:
            clean_name = group_name.replace('(', '').replace(')', '')
            substituents[int(pos)] = clean_name
            parsed_substituents.append(clean_name)
            
        return substituents, parsed_substituents

    def check_rules(substituents, name_order):
        """Checks all constraints and IUPAC rules for a given parsed structure."""
        # Helper to find position by substituent name
        def find_pos(sub_name):
            for pos, name in substituents.items():
                if name == sub_name:
                    return pos
            return None

        # Check for presence of all required groups
        required_groups = ['COOH', 'formyl', 'cyano', 'hydroxy', 'dimethylamino', 'methoxy']
        if not all(group in substituents.values() for group in required_groups):
            return f"FAIL: Not all required substituent groups are present. Found: {list(substituents.values())}"

        pos_cooh = find_pos('COOH')
        pos_formyl = find_pos('formyl')
        pos_cyano = find_pos('cyano')
        pos_hydroxy = find_pos('hydroxy')
        pos_dimethylamino = find_pos('dimethylamino')
        pos_methoxy = find_pos('methoxy')

        # Helper functions for relative positions on a 6-member ring
        def is_meta(p1, p2): return abs(p1 - p2) == 2 or abs(p1 - p2) == 4
        def is_ortho(p1, p2): return abs(p1 - p2) == 1 or abs(p1 - p2) == 5
        def is_para(p1, p2): return abs(p1 - p2) == 3

        # Constraint 1: COOH, formyl, cyano are meta to one another.
        if not (is_meta(pos_cooh, pos_formyl) and is_meta(pos_cooh, pos_cyano) and is_meta(pos_formyl, pos_cyano)):
            return "FAIL: The 'COOH', 'formyl', and 'cyano' groups are not all meta to one another."

        # Constraint 2: hydroxy and dimethylamino are ortho to COOH.
        if not (is_ortho(pos_cooh, pos_hydroxy) and is_ortho(pos_cooh, pos_dimethylamino)):
            return "FAIL: 'hydroxy' and 'dimethylamino' are not both ortho to 'COOH'."

        # Constraint 3: methoxy is para to COOH.
        if not is_para(pos_cooh, pos_methoxy):
            return "FAIL: 'methoxy' is not para to 'COOH'."

        # Constraint 4: methoxy and hydroxy are ortho to cyano.
        if not (is_ortho(pos_cyano, pos_methoxy) and is_ortho(pos_cyano, pos_hydroxy)):
            return "FAIL: 'methoxy' and 'hydroxy' are not both ortho to 'cyano'."
            
        # Constraint 5: IUPAC numbering tie-breaker.
        # Between the ortho groups (dimethylamino and hydroxy), the one first alphabetically
        # ('d' before 'h') gets the lower number.
        if not (pos_dimethylamino < pos_hydroxy):
            return f"FAIL: IUPAC numbering tie-breaker violated. 'dimethylamino' (pos {pos_dimethylamino}) must have a lower number than 'hydroxy' (pos {pos_hydroxy})."

        # Constraint 6: IUPAC alphabetical naming order.
        if name_order != sorted(name_order):
            return f"FAIL: Substituents are not listed in alphabetical order in the name. Order was {name_order}, should be {sorted(name_order)}."

        return "PASS"

    # --- Main execution loop ---
    results = {}
    for key, name in candidates.items():
        substituents, name_order = parse_name(name)
        if substituents is None:
            results[key] = name_order # Store error message
            continue
        
        result_message = check_rules(substituents, name_order)
        results[key] = result_message

    # --- Final Verdict ---
    llm_answer_is_correct = True
    failure_reasons = []

    # Check if the LLM's chosen answer passed all checks
    if results[llm_answer_key] != "PASS":
        llm_answer_is_correct = False
        failure_reasons.append(f"The proposed answer {llm_answer_key} is incorrect. Reason: {results[llm_answer_key]}")

    # Check if any other answer also passed, which would indicate an ambiguous question or faulty check
    for key, result in results.items():
        if key != llm_answer_key and result == "PASS":
            llm_answer_is_correct = False
            failure_reasons.append(f"Another candidate, {key}, also passed all checks, making the provided answer ambiguous or incorrect.")

    if llm_answer_is_correct:
        return "Correct"
    else:
        # Provide a detailed reason for the failure
        return "\n".join(failure_reasons)

# Run the check and print the result
result = check_correctness()
if result == "Correct":
    print("Correct")
else:
    print(f"Incorrect. Reason(s):\n{result}")
