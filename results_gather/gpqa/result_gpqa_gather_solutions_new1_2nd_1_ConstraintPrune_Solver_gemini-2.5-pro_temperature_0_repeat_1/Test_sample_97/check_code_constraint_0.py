import re

def check_chemistry_answer():
    """
    This function checks the correctness of the final answer by applying chemical constraints.

    The reaction is a Ring-Opening Cross-Metathesis (ROCM):
    A + catalyst + 1-propene ---> 1-(prop-1-en-1-yl)-2-vinylcyclopentane

    This imposes several constraints on the starting material 'A':
    1.  Reaction Type: 'A' must be a strained, bicyclic alkene to undergo ring-opening.
    2.  Product Core: 'A' must contain a five-membered ring that remains intact. This means the reactive double bond must be in the *other* ring.
    3.  Product Regiochemistry: The product is 1,2-disubstituted. This means the ring that opens must be fused to two adjacent carbons of the five-membered ring.
    """
    
    # The final answer provided by the LLM to be checked.
    final_answer_text = "<<<A>>>"

    # Extract the letter from the answer format, e.g., 'A' from '<<<A>>>'
    match = re.search(r'<<<([A-D])>>>', final_answer_text)
    if not match:
        return "Invalid answer format. Expected format like <<<A>>>."
    proposed_answer_letter = match.group(1)

    # Define the properties of each candidate based on chemical knowledge.
    candidates = {
        'A': {
            'name': 'bicyclo[3.2.0]hept-6-ene',
            'is_bicyclic': True,
            'has_cyclopentane_core': True,
            'db_in_strained_ring': True,  # Double bond is in the 4-membered ring, not the 5-membered one.
            'fusion_leads_to_1_2_subst': True  # [3.2.0] fusion is at adjacent carbons.
        },
        'B': {
            'name': '1,2-dimethylenecyclopentane',
            'is_bicyclic': False, # It's a monocyclic diene.
            'has_cyclopentane_core': True,
            'db_in_strained_ring': None, # Not applicable.
            'fusion_leads_to_1_2_subst': None # Not applicable.
        },
        'C': {
            'name': '2-methyl-3-methylenebicyclo[2.1.0]pentane',
            'is_bicyclic': True,
            'has_cyclopentane_core': False, # It's a bicyclo[2.1.0]pentane (fused 3 & 4 membered rings).
            'db_in_strained_ring': True,
            'fusion_leads_to_1_2_subst': False # Would not produce a cyclopentane core.
        },
        'D': {
            'name': '2-methylbicyclo[3.1.0]hex-2-ene',
            'is_bicyclic': True,
            'has_cyclopentane_core': True,
            'db_in_strained_ring': False, # Double bond is in the 5-membered ring.
            'fusion_leads_to_1_2_subst': False # Opening the 5-membered ring is incorrect.
        }
    }

    # Determine the logically correct answer by applying the constraints.
    correct_candidate_letter = None
    for letter, properties in candidates.items():
        if (properties['is_bicyclic'] and
            properties['has_cyclopentane_core'] and
            properties['db_in_strained_ring'] and
            properties['fusion_leads_to_1_2_subst']):
            correct_candidate_letter = letter
            break  # Assume only one correct answer.

    # Check if the proposed answer is the correct one.
    if proposed_answer_letter == correct_candidate_letter:
        return "Correct"
    else:
        # Provide a specific reason why the proposed answer is wrong.
        wrong_answer_properties = candidates[proposed_answer_letter]
        
        if not wrong_answer_properties['is_bicyclic']:
            return f"Incorrect. The proposed answer '{wrong_answer_properties['name']}' is not a bicyclic alkene and cannot undergo a Ring-Opening Cross-Metathesis (ROCM) reaction."
            
        if not wrong_answer_properties['has_cyclopentane_core']:
            return f"Incorrect. The proposed answer '{wrong_answer_properties['name']}' does not contain a five-membered ring, so it cannot produce the cyclopentane core of the product."
            
        if not wrong_answer_properties['db_in_strained_ring']:
            return f"Incorrect. In the proposed answer '{wrong_answer_properties['name']}', the double bond is within the five-membered ring. ROCM would open this ring, destroying the required cyclopentane core of the product."
            
        return f"Incorrect. The structure of the proposed answer '{wrong_answer_properties['name']}' does not satisfy all reaction constraints. The correct answer is {correct_candidate_letter} ({candidates[correct_candidate_letter]['name']})."

# Run the check and print the result.
result = check_chemistry_answer()
print(result)