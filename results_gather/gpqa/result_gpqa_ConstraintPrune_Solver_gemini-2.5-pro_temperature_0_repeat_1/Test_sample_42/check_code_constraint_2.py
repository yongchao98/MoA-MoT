import re

def check_iupac_name(name_string: str) -> (bool, str):
    """
    Checks a given IUPAC name against all question constraints and IUPAC rules.
    Returns a tuple (is_fully_correct, reason_string).
    """
    
    # --- 1. Parse the name to get substituent positions ---
    substituents = {1: "COOH"}
    name_order_groups = []
    try:
        substituent_part = name_string.split('benzoic acid')[0]
        pattern = r'(\d+)-(\(?[a-zA-Z]+\)?)'
        matches = re.findall(pattern, substituent_part)
        for pos, group in matches:
            clean_group = group.replace('(', '').replace(')', '')
            substituents[int(pos)] = clean_group
            name_order_groups.append(clean_group)
    except Exception:
        return False, "Failed to parse the name string."

    # Check if all required groups are present in the parsed structure
    required_groups = {"COOH", "formyl", "cyano", "hydroxy", "dimethylamino", "methoxy"}
    if set(substituents.values()) != required_groups:
        return False, f"The name does not contain all the required substituents. Found: {set(substituents.values())}"

    # Helper function to find position of a substituent
    def find_pos(sub_name):
        for pos, name in substituents.items():
            if name == sub_name:
                return pos
        return None

    # --- 2. Check Structural Constraints from the question ---
    pos_cooh, pos_formyl, pos_cyano, pos_hydroxy, pos_dimethylamino, pos_methoxy = (
        1, find_pos("formyl"), find_pos("cyano"), find_pos("hydroxy"),
        find_pos("dimethylamino"), find_pos("methoxy")
    )

    def is_ortho(p1, p2): return abs(p1 - p2) == 1 or abs(p1 - p2) == 5
    def is_meta(p1, p2): return abs(p1 - p2) == 2 or abs(p1 - p2) == 4
    def is_para(p1, p2): return abs(p1 - p2) == 3

    if not (is_meta(pos_cooh, pos_formyl) and is_meta(pos_cooh, pos_cyano) and is_meta(pos_formyl, pos_cyano)):
        return False, "Structural Constraint Failed: COOH, formyl, and cyano are not all meta to one another."
    if not (is_ortho(pos_cooh, pos_hydroxy) and is_ortho(pos_cooh, pos_dimethylamino)):
        return False, "Structural Constraint Failed: Hydroxy and dimethylamino are not both ortho to COOH."
    if not is_para(pos_cooh, pos_methoxy):
        return False, "Structural Constraint Failed: Methoxy is not para to COOH."
    if not (is_ortho(pos_cyano, pos_methoxy) and is_ortho(pos_cyano, pos_hydroxy)):
        return False, "Structural Constraint Failed: Methoxy and hydroxy are not both ortho to the cyano group."

    # --- 3. Check IUPAC Naming Rules ---
    
    # Rule A: Substituents must be listed alphabetically in the name.
    alphabetical_order_groups = sorted(name_order_groups)
    if name_order_groups != alphabetical_order_groups:
        return False, f"IUPAC Rule Failed: Substituents are not listed alphabetically. Order is {name_order_groups}, should be {alphabetical_order_groups}."

    # Rule B: Lowest Locant Tie-breaker for Numbering.
    # When the locant set is fixed (here {2,3,4,5,6}), numbering is chosen to give the lowest locants
    # to the substituents when they are ordered alphabetically.
    alpha_prefixes = sorted([s for s in substituents.values() if s != 'COOH'])
    
    # Get the locants for these prefixes based on the name's numbering
    locants_from_name = [find_pos(p) for p in alpha_prefixes]
    
    # Get the locants for the alternative numbering scheme (for a 6-ring, p -> 8-p)
    locants_alternative = [8 - find_pos(p) for p in alpha_prefixes]
    
    # The locant list from the name must be lexicographically smaller than or equal to the alternative.
    if locants_from_name > locants_alternative:
        return False, f"IUPAC Rule Failed: Incorrect numbering. The locant sequence for alphabetically ordered substituents {locants_from_name} is not the lowest possible. The correct sequence would be {locants_alternative}."

    # If all checks pass, the name is correct.
    return True, "Correct"

# --- Main Execution ---
candidates = {
    "A": "2-(dimethylamino)-3-formyl-4-methoxy-5-cyano-6-hydroxybenzoic acid",
    "B": "5-cyano-2-(dimethylamino)-3-formyl-6-hydroxy-4-methoxybenzoic acid",
    "C": "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid",
    "D": "2-hydroxy-3-cyano-4-methoxy-5-formyl-6-(dimethylamino)benzoic acid"
}

llm_answer_key = "B"
llm_answer_name = candidates[llm_answer_key]

is_correct, reason = check_iupac_name(llm_answer_name)

if is_correct:
    print("Correct")
else:
    print(f"The provided answer '{llm_answer_key}' is incorrect for the following reason:")
    print(f"Reason: {reason}")
    
    print("\n--- Verifying all candidates ---")
    correct_answer_found = False
    for key, name in candidates.items():
        is_cand_correct, reason_cand = check_iupac_name(name)
        if is_cand_correct:
            print(f"Candidate '{key}' is correct.")
            correct_answer_found = True
        else:
            print(f"Candidate '{key}' is incorrect. Reason: {reason_cand}")
    
    if not correct_answer_found:
        print("\nNone of the options are correct.")
