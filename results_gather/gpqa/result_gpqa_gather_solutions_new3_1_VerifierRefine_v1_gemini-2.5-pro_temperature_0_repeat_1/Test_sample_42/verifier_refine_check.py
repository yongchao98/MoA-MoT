import re

def check_iupac_name():
    """
    Checks the correctness of the IUPAC name for the given molecule description.
    """
    # The final answer from the LLM to be checked.
    llm_answer_option = "B"
    options = {
        "A": "5-cyano-2-(dimethylamino)-3-formyl-6-hydroxy-4-methoxybenzoic acid",
        "B": "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid",
        "C": "2-hydroxy-3-cyano-4-methoxy-5-formyl-6-(dimethylamino)benzoic acid",
        "D": "2-(dimethylamino)-3-formyl-4-methoxy-5-cyano-6-hydroxybenzoic acid"
    }
    name_to_check = options.get(llm_answer_option)

    # --- 1. Parse the IUPAC name to build the structure ---
    substituents = {1: 'COOH'}  # Parent is benzoic acid
    name_prefix = name_to_check.replace('benzoic acid', '').strip()
    
    # Regex to find locant-substituent pairs, e.g., "3-cyano", "6-(dimethylamino)"
    matches = re.findall(r'(\d+)-\(?(.*?)\)?(?:-|$)', name_prefix)
    
    if not matches or len(matches) != 5:
        return f"Parsing failed: Expected 5 substituents, but found {len(matches)}."

    parsed_substituents = []
    for locant_str, sub_name in matches:
        parsed_substituents.append((int(locant_str), sub_name))

    name_to_formula = {
        'cyano': 'CN', 'dimethylamino': 'N(CH3)2', 'formyl': 'CHO',
        'hydroxy': 'OH', 'methoxy': 'OCH3'
    }
    
    for locant, sub_name in parsed_substituents:
        if sub_name in name_to_formula:
            substituents[locant] = name_to_formula[sub_name]
        else:
            return f"Parsing failed: Unknown substituent '{sub_name}'."

    if len(substituents) != 6:
        return f"Parsing failed: The name should describe 6 groups (1 parent + 5 substituents)."

    # --- 2. Verify the structure against the question's constraints ---
    def is_ortho(p1, p2): return abs(p1 - p2) == 1 or abs(p1 - p2) == 5
    def is_meta(p1, p2): return abs(p1 - p2) == 2 or abs(p1 - p2) == 4
    def is_para(p1, p2): return abs(p1 - p2) == 3

    pos = {v: k for k, v in substituents.items()}

    # Constraint 1: COOH, CHO, CN are all meta to one another.
    if not (is_meta(pos['COOH'], pos['CHO']) and is_meta(pos['COOH'], pos['CN']) and is_meta(pos['CHO'], pos['CN'])):
        return "Constraint check failed: Carboxylic acid, carbaldehyde, and cyano groups are not all meta to one another."

    # Constraint 2: OH and N(CH3)2 are ortho to COOH.
    if not (is_ortho(pos['COOH'], pos['OH']) and is_ortho(pos['COOH'], pos['N(CH3)2'])):
        return "Constraint check failed: Hydroxyl and dimethylamino groups are not both ortho to the carboxylic acid."

    # Constraint 3: Methoxy is para to COOH.
    if not is_para(pos['COOH'], pos['OCH3']):
        return "Constraint check failed: Methoxy group is not para to the carboxylic acid."

    # Constraint 4: Methoxy and alcohol (OH) are both ortho to the nitrile (CN).
    if not (is_ortho(pos['CN'], pos['OCH3']) and is_ortho(pos['CN'], pos['OH'])):
        return "Constraint check failed: Methoxy and hydroxyl groups are not both ortho to the cyano group."

    # --- 3. Verify IUPAC naming rules for the name itself ---
    # Rule 1: Alphabetical order of substituents in the name string.
    sub_names_in_order = [s[1] for s in parsed_substituents]
    if sub_names_in_order != sorted(sub_names_in_order):
        return f"Naming rule check failed: The substituents are not listed in alphabetical order. Found: {sub_names_in_order}, Expected: {sorted(sub_names_in_order)}."

    # Rule 2: Lowest locant tie-breaker.
    # Generate the alternative numbering (mirror image) and compare locant sets.
    mirror_pos = {group: (loc if loc in [1, 4] else (8 - loc)) for group, loc in pos.items()}
    
    alphabetical_sub_names = sorted(name_to_formula.keys())
    
    original_locant_set = [pos[name_to_formula[name]] for name in alphabetical_sub_names]
    mirror_locant_set = [mirror_pos[name_to_formula[name]] for name in alphabetical_sub_names]

    # The original locant set must be lexicographically smaller than or equal to any alternative.
    if original_locant_set > mirror_locant_set:
        return f"Naming rule check failed: The numbering does not give the lowest locant set based on alphabetical order. Found locant set {original_locant_set}, but the alternative numbering gives a lower set {mirror_locant_set}."

    return "Correct"

# Run the check
result = check_iupac_name()
print(result)