import re

def check_rcm_answer():
    """
    Checks the correctness of the LLM's answer for an RCM synthesis problem.
    """
    # --- Problem Definition ---
    question_options = {
        'A': '5-isopropyl-3,4-dimethylocta-2,6-diene',
        'B': '4-isopropyl-5,6-dimethylocta-1,7-diene',
        'C': '5-isopropyl-3,4-dimethylocta-1,7-diene',
        'D': '5-isopropyl-3,4-dimethylocta-1,6-diene'
    }
    
    target_product_name = "5-isopropyl-3,4-dimethylcyclohex-1-ene"
    # The structure of the target product, based on its IUPAC name.
    target_product_structure = {
        'ring_type': 'cyclohex-1-ene',
        'substituents': {'3': 'methyl', '4': 'methyl', '5': 'isopropyl'}
    }
    
    llm_answer = 'B'

    # --- Helper Functions ---
    def parse_name(name: str) -> dict:
        """Parses the IUPAC name into a structural representation."""
        try:
            parent_match = re.search(r'octa-\d+,\d+-diene', name)
            if not parent_match:
                return {'error': f"Could not find 'octa-x,y-diene' parent in '{name}'"}
            parent = parent_match.group(0)
            
            sub_string = name[:parent_match.start()].strip('-')
            substituents = {}
            parts = sub_string.split('-')
            
            i = 0
            while i < len(parts):
                locants_str = parts[i]
                group_name = parts[i+1]
                base_name = re.sub(r'^(di|tri|tetra)', '', group_name)
                locants = locants_str.split(',')
                for locant in locants:
                    substituents[locant] = base_name
                i += 2
            return {'parent': parent, 'substituents': substituents}
        except Exception:
            return {'error': f"Failed to parse '{name}'"}

    def simulate_rcm(precursor: dict) -> dict:
        """Simulates the RCM reaction and returns the product structure."""
        # Constraint 1: Must be an octa-1,7-diene to form a six-membered ring.
        if precursor.get('parent') != 'octa-1,7-diene':
            return {'error': f"Incorrect precursor type '{precursor.get('parent')}' for 6-membered ring."}

        # Constraint 2: Apply the correct atom mapping from precursor (D) to product (P).
        # This "reversed" mapping is key: P3<-D6, P4<-D5, P5<-D4.
        precursor_to_product_map = {'6': '3', '5': '4', '4': '5', '3': '6'}
        
        product_substituents = {}
        for p_pos, sub_type in precursor.get('substituents', {}).items():
            if p_pos in precursor_to_product_map:
                r_pos = precursor_to_product_map[p_pos]
                product_substituents[r_pos] = sub_type
        
        return {'ring_type': 'cyclohex-1-ene', 'substituents': product_substituents}

    # --- Verification Logic ---
    correct_option = None
    for option_key, option_name in question_options.items():
        precursor_structure = parse_name(option_name)
        if 'error' in precursor_structure:
            continue

        simulated_product = simulate_rcm(precursor_structure)
        if 'error' in simulated_product:
            continue
        
        # Compare the simulated product's substituents with the target's.
        # Sorting items ensures comparison is order-independent.
        if sorted(simulated_product['substituents'].items()) == sorted(target_product_structure['substituents'].items()):
            correct_option = option_key
            break
            
    # --- Final Verdict ---
    if llm_answer == correct_option:
        return "Correct"
    else:
        reason = f"The provided answer is '{llm_answer}', but the correct answer is '{correct_option}'.\n\n"
        reason += f"Analysis:\n"
        reason += f"1. The target product is '{target_product_name}', which has substituents at positions 3 (methyl), 4 (methyl), and 5 (isopropyl).\n"
        reason += "2. To form a 6-membered ring via RCM, the precursor must be an 'octa-1,7-diene'. This eliminates options A and D.\n"
        reason += "3. A forward synthesis simulation shows:\n"
        
        # Analysis of the correct choice
        correct_precursor = parse_name(question_options[correct_option])
        reason += f"   - For option '{correct_option}' ({question_options[correct_option]}):\n"
        reason += f"     - The precursor has substituents at positions {list(correct_precursor['substituents'].keys())}.\n"
        reason += f"     - Applying the correct RCM mapping (P3<-D6, P4<-D5, P5<-D4) correctly places the substituents at positions 3, 4, and 5 on the product ring, matching the target.\n"
        
        # Analysis of the LLM's choice (if it was wrong and a valid precursor type)
        if llm_answer in ['B', 'C']:
            llm_precursor = parse_name(question_options[llm_answer])
            llm_product = simulate_rcm(llm_precursor)
            reason += f"   - For option '{llm_answer}' ({question_options[llm_answer]}):\n"
            reason += f"     - The precursor has substituents at positions {list(llm_precursor['substituents'].keys())}.\n"
            reason += f"     - This results in a product with substituents at positions {list(llm_product['substituents'].keys())}, which does not match the target.\n"
            
        return reason

# Run the check and print the result.
result = check_rcm_answer()
print(result)