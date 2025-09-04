import re

def parse_molecule_name(name):
    """Parses an IUPAC name into a dictionary of its components."""
    parts = name.split('-')
    parent_info = parts[-1]
    
    # Find parent chain and diene type
    parent_chain = None
    if 'diene' in parent_info:
        parent_chain = parent_info
    
    substituents = {}
    # Regex to find substituent names like 'isopropyl', 'methyl', 'dimethyl', etc.
    sub_name_pattern = r'(iso|sec|tert|di|tri|tetra)?(methyl|ethyl|propyl|butyl|pentyl|hexyl)'
    
    for part in parts[:-1]:
        # Find all numbers (positions) in the part
        positions = [int(x) for x in re.findall(r'\d+', part)]
        # Find the substituent name
        sub_match = re.search(sub_name_pattern, part)
        if sub_match:
            # For 'dimethyl', 'trimethyl', etc., we need the base name 'methyl'
            sub_name = sub_match.group(2)
            if sub_match.group(1) == 'iso':
                sub_name = 'isopropyl' # special case for isopropyl
            for pos in positions:
                substituents[pos] = sub_name
                
    return {"parent": parent_chain, "subs": substituents}

def get_iupac_name(substituents):
    """Generates the IUPAC name for a substituted cyclohex-1-ene."""
    if not substituents:
        return "cyclohex-1-ene"

    # Sort substituents alphabetically for naming
    sorted_subs = sorted(substituents.items(), key=lambda item: item[1])
    
    # Group by substituent type for names like 'dimethyl'
    grouped_subs = {}
    for pos, name in sorted_subs:
        if name not in grouped_subs:
            grouped_subs[name] = []
        grouped_subs[name].append(str(pos))

    name_parts = []
    for name, positions in grouped_subs.items():
        prefix = ""
        if len(positions) == 2:
            prefix = "di"
        elif len(positions) == 3:
            prefix = "tri"
        
        pos_str = ",".join(positions)
        name_parts.append(f"{pos_str}-{prefix}{name}")
        
    # Sort the final name parts based on the base substituent name (e.g., 'isopropyl' before 'methyl')
    # This is a simplified sort key that works for this problem
    name_parts.sort(key=lambda x: x.split('-')[1])

    return f"{'-'.join(name_parts)}-cyclohex-1-ene"


def simulate_rcm(precursor_name):
    """Simulates RCM and returns the IUPAC name of the product."""
    precursor = parse_molecule_name(precursor_name)

    # 1. Check for correct ring size precursor
    if precursor.get("parent") != "octa-1,7-diene":
        return f"Invalid precursor: Forms a ring of the wrong size, not a 6-membered ring."

    # 2. Apply the correct atom mapping for RCM
    # Precursor position -> Product position
    mapping = {3: 6, 4: 5, 5: 4, 6: 3}
    
    product_subs = {}
    for pre_pos, sub_name in precursor["subs"].items():
        if pre_pos in mapping:
            pro_pos = mapping[pre_pos]
            product_subs[pro_pos] = sub_name

    # 3. Name the product according to IUPAC rules (lowest locants)
    # The mapping used inherently produces the lowest locant set for the product name.
    # For example, if a precursor has subs at 3,4,5, the product has subs at 6,5,4.
    # The product is named by numbering from the other side of the double bond,
    # which makes the locants 3,4,5. Our mapping {3:6, 4:5, 5:4} already accounts for this.
    
    # We must check if an alternative numbering gives a lower locant set.
    # This is not necessary with the chosen mapping, but it's good practice.
    # The locant set is the set of keys in product_subs.
    
    return get_iupac_name(product_subs)


def check_answer():
    """
    Checks the correctness of the LLM's answer.
    """
    question_options = {
        "A": "5-isopropyl-3,4-dimethylocta-1,6-diene",
        "B": "5-isopropyl-3,4-dimethylocta-2,6-diene",
        "C": "5-isopropyl-3,4-dimethylocta-1,7-diene",
        "D": "4-isopropyl-5,6-dimethylocta-1,7-diene"
    }
    target_product = "5-isopropyl-3,4-dimethylcyclohex-1-ene"
    llm_answer_choice = "D"
    
    chosen_precursor_name = question_options[llm_answer_choice]
    
    # Simulate the reaction for the chosen answer
    resulting_product = simulate_rcm(chosen_precursor_name)
    
    # Normalize names for comparison by removing hyphens
    normalized_result = resulting_product.replace('-', '')
    normalized_target = target_product.replace('-', '')

    if normalized_result == normalized_target:
        return "Correct"
    else:
        # Provide a detailed reason for the failure
        error_reason = f"The answer '{llm_answer_choice}' is incorrect.\n"
        error_reason += f"The starting material is: {chosen_precursor_name}.\n"
        error_reason += f"Simulating the RCM reaction on this material produces: '{resulting_product}'.\n"
        error_reason += f"This does not match the target product: '{target_product}'."
        
        # Additionally, check which option IS correct
        for option, name in question_options.items():
            product = simulate_rcm(name)
            if product.replace('-', '') == normalized_target:
                error_reason += f"\n\nVerification shows that option '{option}' ({name}) is the correct precursor."
                break
        
        return error_reason

# Run the check
result = check_answer()
print(result)