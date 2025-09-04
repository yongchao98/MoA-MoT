import re

def check_iupac_name():
    """
    Checks the correctness of the IUPAC name for the described molecule.
    """
    llm_provided_answer = "D"
    options = {
        "A": "5-cyano-2-(dimethylamino)-3-formyl-6-hydroxy-4-methoxybenzoic acid",
        "B": "2-hydroxy-3-cyano-4-methoxy-5-formyl-6-(dimethylamino)benzoic acid",
        "C": "2-(dimethylamino)-3-formyl-4-methoxy-5-cyano-6-hydroxybenzoic acid",
        "D": "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid"
    }

    # --- Step 1: Define helper functions and possible structures ---

    def is_ortho(p1, p2): return abs(p1 - p2) == 1 or {p1, p2} == {1, 6}
    def is_meta(p1, p2): return abs(p1 - p2) == 2 or {p1, p2} == {1, 5} or {p1, p2} == {2, 6}
    def is_para(p1, p2): return abs(p1 - p2) == 3

    # There are two potential structures that fit the geometric description.
    # Structure 1 (counter-clockwise numbering from COOH)
    struct1 = {1: 'COOH', 2: 'N(CH3)2', 3: 'CHO', 4: 'OCH3', 5: 'CN', 6: 'OH'}
    # Structure 2 (clockwise numbering from COOH)
    struct2 = {1: 'COOH', 2: 'OH', 3: 'CN', 4: 'OCH3', 5: 'CHO', 6: 'N(CH3)2'}

    # --- Step 2: Verify both structures meet all constraints ---

    def verify_constraints(structure, struct_name):
        pos = {v: k for k, v in structure.items()}
        # C1: COOH, CHO, CN are meta
        if not (is_meta(pos['COOH'], pos['CHO']) and is_meta(pos['COOH'], pos['CN']) and is_meta(pos['CHO'], pos['CN'])):
            return f"{struct_name} fails: COOH, CHO, CN not all meta."
        # C2: OH, N(CH3)2 are ortho to COOH
        if not (is_ortho(pos['COOH'], pos['OH']) and is_ortho(pos['COOH'], pos['N(CH3)2'])):
            return f"{struct_name} fails: OH/N(CH3)2 not ortho to COOH."
        # C3: OCH3 is para to COOH
        if not is_para(pos['COOH'], pos['OCH3']):
            return f"{struct_name} fails: OCH3 not para to COOH."
        # C4: OCH3, OH are ortho to CN
        if not (is_ortho(pos['CN'], pos['OCH3']) and is_ortho(pos['CN'], pos['OH'])):
            return f"{struct_name} fails: OCH3/OH not ortho to CN."
        return "Valid"

    if verify_constraints(struct1, "Structure 1") != "Valid": return verify_constraints(struct1, "Structure 1")
    if verify_constraints(struct2, "Structure 2") != "Valid": return verify_constraints(struct2, "Structure 2")
    # Both structures are geometrically valid based on the description.

    # --- Step 3: Apply IUPAC numbering tie-breaker rule ---

    substituent_prefixes = {
        'CN': 'cyano', 'N(CH3)2': 'dimethylamino', 'CHO': 'formyl',
        'OH': 'hydroxy', 'OCH3': 'methoxy'
    }
    alphabetical_order = sorted(substituent_prefixes.values())
    first_alpha_sub_prefix = alphabetical_order[0] # 'cyano'
    first_alpha_sub_code = 'CN'

    # Find the locant for the first alphabetical substituent in each structure
    locant_in_struct1 = {v: k for k, v in struct1.items()}[first_alpha_sub_code] # is 5
    locant_in_struct2 = {v: k for k, v in struct2.items()}[first_alpha_sub_code] # is 3

    # The correct numbering gives the lowest locant at the first point of difference.
    # Since 3 < 5, Structure 2 is the correctly numbered one.
    if locant_in_struct2 < locant_in_struct1:
        correct_structure = struct2
    else:
        correct_structure = struct1

    # --- Step 4: Construct the correct name from the correct structure ---

    correct_locants = {v: k for k, v in correct_structure.items()}
    name_parts = []
    for prefix in alphabetical_order:
        sub_code = [code for code, p in substituent_prefixes.items() if p == prefix][0]
        locant = correct_locants[sub_code]
        # Handle parentheses for complex substituents
        part = f"({prefix})" if prefix == "dimethylamino" else prefix
        name_parts.append(f"{locant}-{part}")
    
    generated_correct_name = "-".join(name_parts) + "benzoic acid"

    # --- Step 5: Verify the LLM's answer ---
    
    llm_choice_name = options.get(llm_provided_answer)
    if not llm_choice_name:
        return f"LLM answer '{llm_provided_answer}' is not a valid option."

    # Check 5a: Does the name from the answer correspond to the correct structure?
    # We can check this by comparing the generated name with the answer's name.
    if llm_choice_name != generated_correct_name:
        return (f"Incorrect. The IUPAC name in answer {llm_provided_answer} is wrong. "
                f"Based on numbering rules (lowest locant for 'cyano'), the correct structure is {correct_structure}. "
                f"The correct alphabetical name for this structure is '{generated_correct_name}'.")

    # Check 5b: Is the name in the answer alphabetized correctly?
    # (This is implicitly checked above, but we can be explicit)
    sub_parts_in_answer = re.findall(r'\d+-(\(?[a-zA-Z]+\)?)', llm_choice_name)
    alpha_key = lambda s: s.replace('(', '').replace(')', '')
    if [alpha_key(s) for s in sub_parts_in_answer] != sorted([alpha_key(s) for s in sub_parts_in_answer]):
        return f"Incorrect. The substituents in answer {llm_provided_answer} are not in alphabetical order."

    return "Correct"

# Execute the check and print the result
result = check_iupac_name()
if result == "Correct":
    print("Correct")
else:
    print(f"Incorrect. Reason: {result}")
