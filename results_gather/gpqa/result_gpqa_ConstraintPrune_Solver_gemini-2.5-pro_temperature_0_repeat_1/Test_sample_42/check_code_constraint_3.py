import re

def get_structure_from_name(name):
    """Parses an IUPAC name string and returns a dictionary representing the molecule's structure."""
    # The parent is benzoic acid, so COOH is at position 1.
    structure = {1: "COOH"}
    # Extract the part of the name with the substituents
    substituent_part = name.split('benzoic acid')[0].strip()
    # Regex to find number-group pairs, handling parentheses for complex groups
    pattern = r'(\d+)-(\(?[a-zA-Z]+\)?)'
    matches = re.findall(pattern, substituent_part)
    for pos, group in matches:
        # Clean up the group name by removing parentheses
        group_name = group.replace('(', '').replace(')', '')
        # Map known prefixes to standard group names
        group_map = {
            "formyl": "CHO",
            "cyano": "CN",
            "hydroxy": "OH",
            "dimethylamino": "N(CH3)2",
            "methoxy": "OCH3"
        }
        structure[int(pos)] = group_map.get(group_name, group_name)
    return structure

def check_correctness(llm_answer_key, all_options):
    """
    Checks the correctness of the provided answer by verifying all constraints and IUPAC rules.
    """
    llm_answer_name = all_options.get(llm_answer_key)
    if not llm_answer_name:
        return f"Invalid answer key '{llm_answer_key}'. It was not found in the options."

    # --- Define Helper Functions for Positional Relationships ---
    def is_ortho(p1, p2): return abs(p1 - p2) % 6 == 1 or abs(p1 - p2) % 6 == 5
    def is_meta(p1, p2): return abs(p1 - p2) % 6 == 2 or abs(p1 - p2) % 6 == 4
    def is_para(p1, p2): return abs(p1 - p2) % 6 == 3

    # --- Define the Correct Locant Set based on IUPAC rules ---
    # Alphabetical order: cyano, dimethylamino, formyl, hydroxy, methoxy
    # The lowest possible locant set for these substituents is [3, 6, 5, 2, 4]
    correct_locant_set = [3, 6, 5, 2, 4]

    # --- Check each option ---
    for key, name in all_options.items():
        is_correct_choice = (key == llm_answer_key)
        
        # 1. Parse the name to get the structure
        try:
            structure = get_structure_from_name(name)
            # Helper to find position of a group
            def find_pos(group_name):
                for pos, name_in_struct in structure.items():
                    if name_in_struct == group_name:
                        return pos
                return None
            
            pos_cooh, pos_cho, pos_cn, pos_oh, pos_nme2, pos_ome = (
                find_pos('COOH'), find_pos('CHO'), find_pos('CN'),
                find_pos('OH'), find_pos('N(CH3)2'), find_pos('OCH3')
            )
            
            if None in [pos_cooh, pos_cho, pos_cn, pos_oh, pos_nme2, pos_ome]:
                 if is_correct_choice: return f"The proposed answer '{key}' is incorrect. It is missing one or more required substituent groups."
                 continue

        except Exception as e:
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. Failed to parse the name: {e}"
            continue

        # 2. Check Structural Constraints
        if not (is_meta(pos_cooh, pos_cho) and is_meta(pos_cooh, pos_cn) and is_meta(pos_cho, pos_cn)):
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. The COOH, formyl, and cyano groups are not all meta to one another."
            continue
        if not (is_ortho(pos_cooh, pos_oh) and is_ortho(pos_cooh, pos_nme2)):
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. The hydroxyl and dimethylamino groups are not both ortho to the carboxylic acid."
            continue
        if not is_para(pos_cooh, pos_ome):
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. The methoxy group is not para to the carboxylic acid."
            continue
        if not (is_ortho(pos_cn, pos_ome) and is_ortho(pos_cn, pos_oh)):
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. The methoxy and hydroxyl groups are not both ortho to the cyano group."
            continue

        # 3. Check IUPAC Naming Rules
        # Rule 3a: Alphabetical order of substituents in the name string
        substituent_part = name.split('benzoic acid')[0].strip()
        matches = re.findall(r'\d+-(\(?[a-zA-Z]+\)?)', substituent_part)
        groups_in_name = [g.replace('(', '').replace(')', '') for g in matches]
        sorted_groups = sorted(groups_in_name)
        if groups_in_name != sorted_groups:
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. The substituents are not listed in alphabetical order. Found '{groups_in_name}', expected '{sorted_groups}'."
            continue

        # Rule 3b: Lowest Locant Set
        substituent_names_alpha_order = ['cyano', 'dimethylamino', 'formyl', 'hydroxy', 'methoxy']
        group_map_reverse = {"CN": "cyano", "N(CH3)2": "dimethylamino", "CHO": "formyl", "OH": "hydroxy", "OCH3": "methoxy"}
        
        candidate_locant_set = []
        for sub_name in substituent_names_alpha_order:
            for pos, group_code in structure.items():
                if group_map_reverse.get(group_code) == sub_name:
                    candidate_locant_set.append(pos)
                    break
        
        if candidate_locant_set != correct_locant_set:
            if is_correct_choice: return f"The proposed answer '{key}' is incorrect. It does not use the lowest possible locant set. Found {candidate_locant_set}, but the correct set is {correct_locant_set}."
            continue
            
        # If all checks pass for a given option, it's a valid name.
        # If this is the LLM's answer, we can break and confirm it's correct.
        if is_correct_choice:
            return "Correct"

    # This part is reached if the provided answer failed a check, but we continued the loop.
    return f"The proposed answer '{llm_answer_key}' is incorrect based on the checks above."


# --- Execution ---
# The question's options
options = {
    "A": "2-(dimethylamino)-3-formyl-4-methoxy-5-cyano-6-hydroxybenzoic acid",
    "B": "5-cyano-2-(dimethylamino)-3-formyl-6-hydroxy-4-methoxybenzoic acid",
    "C": "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid",
    "D": "2-hydroxy-3-cyano-4-methoxy-5-formyl-6-(dimethylamino)benzoic acid"
}

# The answer provided by the other LLM
llm_provided_answer = "C"

# Run the checker
result = check_correctness(llm_provided_answer, options)
print(result)