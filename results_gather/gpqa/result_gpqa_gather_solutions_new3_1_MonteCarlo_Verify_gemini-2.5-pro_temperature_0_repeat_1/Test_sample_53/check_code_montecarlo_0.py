import re

def get_substituents(name):
    """Parses an IUPAC name to find substituents and their positions."""
    substituents = {}
    # Regex to find substituent-position pairs like '5-isopropyl' or '3,4-dimethyl'
    # It handles simple alkyl groups.
    matches = re.findall(r'(\d+(?:,\d+)*)-([a-zA-Z]+)', name)
    for positions_str, sub_name in matches:
        positions = [int(p) for p in positions_str.split(',')]
        
        # Handle prefixes like 'di', 'tri'
        if sub_name.startswith('di'):
            base_name = sub_name[2:]
        elif sub_name.startswith('tri'):
            base_name = sub_name[3:]
        else:
            base_name = sub_name
            
        # Map IUPAC names to a simpler representation
        sub_map = {'methyl': 'Me', 'isopropyl': 'iPr', 'ethyl': 'Et', 'propyl': 'Pr'}
        simple_name = sub_map.get(base_name, base_name)

        for pos in positions:
            substituents[pos] = simple_name
            
    return substituents

def predict_rcm_product_name(precursor_name):
    """
    Predicts the IUPAC name of the product from an RCM reaction of a given diene.
    This function is simplified for octa-1,7-dienes forming cyclohexenes.
    """
    # 1. Check if the precursor is suitable for forming a 6-membered ring
    if 'octa-1,7-diene' not in precursor_name:
        return "Incorrect precursor type: Does not form a six-membered ring."

    # 2. Get substituents from the precursor name
    precursor_subs = get_substituents(precursor_name)
    
    # 3. The RCM of an octa-1,7-diene forms a ring from C2-C7.
    # We need to determine the IUPAC name of the resulting substituted cyclohexene.
    # The new double bond is between original C2 and C7.
    
    # Let's check the two possible numbering directions for the product ring.
    # Direction 1: P1=C7, P2=C2. Ring atoms map: P3->C3, P4->C4, P5->C5, P6->C6
    subs1 = {p: precursor_subs.get(p) for p in range(3, 7) if precursor_subs.get(p)}
    locants1 = sorted(subs1.keys())

    # Direction 2: P1=C2, P2=C7. Ring atoms map: P3->C6, P4->C5, P5->C4, P6->C3
    subs2 = {}
    for p_new in range(3, 7):
        p_old = 9 - p_new # Mapping: 3->6, 4->5, 5->4, 6->3
        if p_old in precursor_subs:
            subs2[p_new] = precursor_subs[p_old]
    locants2 = sorted(subs2.keys())

    # 4. Choose the direction that gives the lowest locant set
    final_subs = {}
    if tuple(locants1) < tuple(locants2):
        final_subs = subs1
    elif tuple(locants2) < tuple(locants1):
        final_subs = subs2
    else: # Tie-break with alphabetization (simplified check)
        # This level of detail is often where errors occur, but for this problem,
        # the locant sets are distinct. We'll assume the first valid one.
        final_subs = subs1 if locants1 else subs2

    if not final_subs:
        return "cyclohex-1-ene"

    # 5. Construct the final product name
    # Invert the dictionary to group by substituent type
    sub_groups = {}
    for pos, name in final_subs.items():
        if name not in sub_groups:
            sub_groups[name] = []
        sub_groups[name].append(str(pos))
    
    prefix_map = {1: '', 2: 'di', 3: 'tri', 4: 'tetra'}
    name_map = {'Me': 'methyl', 'iPr': 'isopropyl'}

    parts = []
    # Sort alphabetically by substituent name for correct IUPAC order
    for sub_simple_name in sorted(sub_groups.keys(), key=lambda s: name_map.get(s,s)):
        positions = sorted(sub_groups[sub_simple_name])
        prefix = prefix_map[len(positions)]
        full_sub_name = name_map.get(sub_simple_name, sub_simple_name)
        parts.append(f"{','.join(positions)}-{prefix}{full_sub_name}")

    return f"{'-'.join(parts)}cyclohex-1-ene".replace('--', '-')


def check_answer():
    """
    Checks the correctness of the LLM's answer for the RCM question.
    """
    question = {
        "target_product": "5-isopropyl-3,4-dimethylcyclohex-1-ene",
        "options": {
            "A": "5-isopropyl-3,4-dimethylocta-2,6-diene",
            "B": "5-isopropyl-3,4-dimethylocta-1,6-diene",
            "C": "5-isopropyl-3,4-dimethylocta-1,7-diene",
            "D": "4-isopropyl-5,6-dimethylocta-1,7-diene"
        },
        "llm_answer": "C"
    }

    chosen_option_name = question["options"][question["llm_answer"]]
    
    # Predict the product from the chosen starting material
    predicted_product = predict_rcm_product_name(chosen_option_name)

    # The predicted name might have substituents in a different order but be chemically identical.
    # A robust check would parse both names into substituent dictionaries and compare them.
    target_subs = get_substituents(question["target_product"])
    predicted_subs = get_substituents(predicted_product)

    if target_subs == predicted_subs:
        # The chosen answer is correct. Now, let's explain why the others are wrong.
        return "Correct"
    else:
        # Analyze why the chosen answer was wrong
        if "Incorrect precursor type" in predicted_product:
            return (f"Incorrect. The chosen starting material, {chosen_option_name}, "
                    f"is of an incorrect type. It would not form a six-membered ring via RCM.")
        else:
            return (f"Incorrect. The chosen starting material, {chosen_option_name}, "
                    f"would produce '{predicted_product}', which is not the target product "
                    f"'{question['target_product']}'. The substituents are in the wrong positions.")

# Run the check
result = check_answer()
print(result)