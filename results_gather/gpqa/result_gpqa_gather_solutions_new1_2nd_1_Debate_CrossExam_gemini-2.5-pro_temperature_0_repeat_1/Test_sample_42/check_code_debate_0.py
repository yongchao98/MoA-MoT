import re

def parse_name(name):
    """Parses an IUPAC name to extract substituents and their locants."""
    base_name = name.replace("benzoic acid", "").strip()
    # Regex to find parts like "3-cyano" or "6-(dimethylamino)"
    parts = re.findall(r'(\d+-\(?[a-zA-Z]+\)?)', base_name)
    substituents = {}
    for part in parts:
        locant, sub_name = part.split('-', 1)
        sub_name = sub_name.strip('()')
        substituents[sub_name] = int(locant)
    return substituents

def check_structural_constraints(substituents):
    """Checks if a parsed structure matches all constraints from the question."""
    # Constraint: para to COOH (C1) is methoxy -> methoxy at C4
    if substituents.get("methoxy") != 4:
        return f"Constraint violated: Methoxy group is not para to the carboxylic acid (should be at C4, but is at {substituents.get('methoxy')})."

    # Constraint: ortho to COOH are hydroxyl and dimethylamino -> at C2 and C6
    ortho_groups = {"hydroxy", "dimethylamino"}
    ortho_positions = {substituents.get("hydroxy"), substituents.get("dimethylamino")}
    if ortho_positions != {2, 6}:
        return f"Constraint violated: Hydroxyl and dimethylamino groups are not ortho to the carboxylic acid (should be at C2 and C6, but are at {ortho_positions})."

    # Constraint: meta to COOH are formyl and cyano -> at C3 and C5
    meta_groups = {"formyl", "cyano"}
    meta_positions = {substituents.get("formyl"), substituents.get("cyano")}
    if meta_positions != {3, 5}:
        return f"Constraint violated: Formyl and cyano groups are not meta to the carboxylic acid (should be at C3 and C5, but are at {meta_positions})."

    # Constraint: methoxy and alcohol (hydroxy) are ortho to nitrile (cyano)
    cyano_pos = substituents.get("cyano")
    hydroxy_pos = substituents.get("hydroxy")
    # Ortho positions on a 6-membered ring, 1-indexed
    ortho_to_cyano = {((cyano_pos - 2 + 6) % 6) + 1, (cyano_pos % 6) + 1}
    
    if 4 not in ortho_to_cyano:
        return f"Constraint violated: Methoxy group (at C4) is not ortho to the cyano group (at C{cyano_pos})."
    if hydroxy_pos not in ortho_to_cyano:
        return f"Constraint violated: Hydroxyl group (at C{hydroxy_pos}) is not ortho to the cyano group (at C{cyano_pos})."

    return "Correct"

def check_alphabetical_order(name):
    """Checks if substituents are listed alphabetically in the name."""
    base_name = name.replace("benzoic acid", "").strip()
    # Extract just the names of the substituents
    sub_names = [re.sub(r'^\d+-', '', part).strip('()') for part in base_name.split('-') if part]
    if sub_names != sorted(sub_names):
        return f"IUPAC rule violated: Substituents are not in alphabetical order. Expected {sorted(sub_names)}, but got {sub_names}."
    return "Correct"

def check_lowest_locant_rule(substituents):
    """Checks if the numbering follows the lowest locant tie-breaker rule."""
    # Generate the alternative (mirror-image) structure
    pos_to_group = {v: k for k, v in substituents.items()}
    alt_substituents = {
        "methoxy": 4,
        pos_to_group[2]: 6,
        pos_to_group[6]: 2,
        pos_to_group[3]: 5,
        pos_to_group[5]: 3,
    }

    # Check if the alternative structure is also valid
    if not check_structural_constraints(alt_substituents) == "Correct":
        # If the alternative is not valid, the proposed one is correct by default.
        return "Correct"

    # If both structures are valid, apply the tie-breaker rule
    alphabetical_sub_order = sorted(substituents.keys())
    
    locants_proposed = [substituents[sub] for sub in alphabetical_sub_order]
    locants_alternative = [alt_substituents[sub] for sub in alphabetical_sub_order]

    if locants_proposed > locants_alternative:
        return f"IUPAC rule violated: Lowest locant rule. The numbering scheme resulting in locant set {locants_alternative} is preferred over {locants_proposed}."
    
    return "Correct"

def check_answer():
    """Main function to check the correctness of the final answer."""
    # The final answer provided by the LLM
    final_answer_name = "3-cyano-6-(dimethylamino)-5-formyl-2-hydroxy-4-methoxybenzoic acid"
    
    # 1. Parse the name to get the structure
    try:
        substituents = parse_name(final_answer_name)
        if len(substituents) != 5:
             return "Parsing Error: Could not correctly parse all 5 substituents from the name."
    except Exception as e:
        return f"Parsing Error: {e}"

    # 2. Check if the structure matches all constraints
    structure_check = check_structural_constraints(substituents)
    if structure_check != "Correct":
        return f"Incorrect. Reason: {structure_check}"

    # 3. Check if the name follows the alphabetical order rule
    alpha_check = check_alphabetical_order(final_answer_name)
    if alpha_check != "Correct":
        return f"Incorrect. Reason: {alpha_check}"

    # 4. Check if the name follows the lowest locant tie-breaker rule
    locant_check = check_lowest_locant_rule(substituents)
    if locant_check != "Correct":
        return f"Incorrect. Reason: {locant_check}"

    # If all checks pass
    return "Correct"

# Run the check
result = check_answer()
print(result)