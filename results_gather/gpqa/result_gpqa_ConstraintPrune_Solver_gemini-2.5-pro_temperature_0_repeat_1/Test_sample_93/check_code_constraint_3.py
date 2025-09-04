import re

def get_precursor_features(name: str) -> dict:
    """
    Parses an IUPAC name to extract features relevant to a Diels-Alder reaction.
    Returns a dictionary with features or an error message.
    """
    # Constraint 1: Reaction Type (Intermolecular vs. Intramolecular)
    if ' and ' in name:
        return {"type": "intermolecular"}

    features = {"type": "intramolecular"}

    # Extract chain length
    chain_map = {'tetradeca': 14, 'trideca': 13, 'dodeca': 12}
    for chain_name, length in chain_map.items():
        if chain_name in name:
            features['chain_length'] = length
            break
    if 'chain_length' not in features:
        return {"error": "Could not determine chain length."}

    # Extract double bond positions to find diene and dienophile
    pos_match = re.search(r'(\d+),(\d+),(\d+)-trien', name)
    if not pos_match:
        return {"error": "Could not determine triene positions."}
    
    positions = sorted([int(p) for p in pos_match.groups()])
    
    # A conjugated diene has positions n, n+2
    if positions[1] == positions[0] + 2:
        features['diene_pos'] = (positions[0], positions[1] + 1)
        features['dienophile_pos'] = (positions[2], positions[2] + 1)
    elif positions[2] == positions[1] + 2:
        features['diene_pos'] = (positions[1], positions[2] + 1)
        features['dienophile_pos'] = (positions[0], positions[0] + 1)
    else:
        return {"error": "No conjugated diene found."}

    # Constraint 2: Linker Length
    # The linker is the set of carbons between the diene and dienophile
    linker_start = min(features['diene_pos'][1], features['dienophile_pos'][1]) + 1
    linker_end = max(features['diene_pos'][0], features['dienophile_pos'][0]) - 1
    features['linker_length'] = linker_end - linker_start + 1

    # Constraint 3: Substituents
    # The ester is at C1. We check the alkyl group at the end of the chain.
    last_reacting_carbon = max(features['diene_pos'][1], features['dienophile_pos'][1])
    substituent_carbons = features['chain_length'] - last_reacting_carbon
    substituent_map = {1: 'methyl', 2: 'ethyl', 3: 'propyl'}
    features['substituent'] = substituent_map.get(substituent_carbons, 'other')
    
    return features

def check_correctness():
    """
    Checks the correctness of the LLM's answer by applying chemical constraints.
    """
    llm_answer = "A"
    candidates = {
        "A": "methyl (2E,4E,10Z)-tetradeca-2,4,10-trienoate",
        "B": "1-vinylcyclohex-1-ene and methyl hex-2-ynoate",
        "C": "Cyclohexene and methyl 2,3-dimethylenehexanoate",
        "D": "methyl (2E,8E,10E)-tetradeca-2,8,10-trienoate"
    }

    # Define target properties based on the product name
    target_properties = {
        "type": "intramolecular",
        "linker_length": 4,
        "substituent": "propyl",
        "diene_pos": (2, 5) # Diene at C2-C5 gives double bond at C3-C4
    }

    for key, name in candidates.items():
        features = get_precursor_features(name)

        # Check against constraints
        if features.get("type") != target_properties["type"]:
            if key == llm_answer:
                return f"Incorrect. The answer {key} is an intermolecular reaction, but an intramolecular reaction is required."
            continue # Correctly fails this candidate

        if features.get("linker_length") != target_properties["linker_length"]:
            return f"Incorrect. Candidate {key} has a linker of length {features.get('linker_length')}, but a 4-atom linker is required."

        if features.get("substituent") != target_properties["substituent"]:
            return f"Incorrect. Candidate {key} provides an '{features.get('substituent')}' group, but the target requires a '{target_properties['substituent']}' group."

        if features.get("diene_pos") != target_properties["diene_pos"]:
            if key == llm_answer:
                 return f"Incorrect. The answer {key} has its diene at C{features.get('diene_pos')[0]}-C{features.get('diene_pos')[1]}, which will not produce the required C3-C4 double bond in the product."
            continue # Correctly fails this candidate
        
        # If a candidate passes all checks, it must be the answer
        if key == llm_answer:
            return "Correct"
        else:
            # This case would mean a different candidate was the correct one
            return f"Incorrect. The correct answer should be {key}, not {llm_answer}."

    # This is reached if the correct answer 'A' failed a check, which is handled above.
    # Or if a non-answer candidate passed, also handled above.
    # This part of the code should ideally not be reached if logic is sound.
    return "An unexpected error occurred in the checking logic."

# Execute the check
result = check_correctness()
print(result)