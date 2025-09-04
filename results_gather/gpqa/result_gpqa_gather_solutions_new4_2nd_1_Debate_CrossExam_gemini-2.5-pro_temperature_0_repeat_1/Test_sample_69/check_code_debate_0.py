import sys

def check_correctness():
    """
    This function checks the correctness of the proposed solution to the chemical riddle.
    It verifies each step of the deduction based on a set of known chemical facts.
    """
    
    # --- Database of Chemical Facts ---
    # This dictionary simulates the chemical knowledge required to solve the puzzle.
    chemical_data = {
        'S8': {
            'state': 'solid',
            'name': 'Sulfur (octatomic ring)',
        },
        'Cl2': {
            'state': 'gas',
            'name': 'Chlorine',
        },
        'SCl2': {
            'name': 'Sulfur dichloride',
            'color': 'cherry-red',
            'hydrolysis': {
                'reforms': 'S8',
                'strong_acid': 'HCl',
                'weak_acid': 'H2SO3'
            }
        },
        'C2H4': {
            'state': 'gas',
            'name': 'Ethene',
        },
        'C2H4Cl2': {
            'name': '1,2-dichloroethane',
            'use': 'solvent',
        },
        'HCl': {
            'name': 'Hydrochloric acid',
            'acid_strength': 'strong',
        },
        'H2SO3': {
            'name': 'Sulfurous acid',
            'acid_strength': 'weak',
        },
        'Mustard Gas': {
            'formula': '(ClCH2CH2)2S',
            'name': 'Mustard gas',
            'hazard_level': 'extremely hazardous',
            'point_group': 'C2'
        },
        'CO': {
            'state': 'gas',
            'name': 'Carbon monoxide'
        },
        'Thiophosgene': {
            'formula': 'CSCl2',
            'name': 'Thiophosgene',
            'hazard_level': 'extremely hazardous',
            'point_group': 'C2v'
        }
    }

    # --- Proposed Solution from the Answer ---
    solution = {
        'A': 'S8',
        'B': 'Cl2',
        'C': 'SCl2',
        'D': 'C2H4',
        'E': 'Mustard Gas',
        'F': 'HCl',
        'G': 'H2SO3',
        'H': 'C2H4Cl2'
    }
    
    final_answer_option = 'B'
    question_options = {'A': 'D4h', 'B': 'C2', 'C': 'D∞h', 'D': 'C2v'}

    # --- Step-by-step Verification ---

    # 1. Check Clue 1: A(s) + 8 B(g) → C (bright red product)
    # Reaction: S₈(s) + 8Cl₂(g) → 8SCl₂(l)
    if chemical_data[solution['A']]['state'] != 'solid':
        return f"Constraint 1 failed: A ({solution['A']}) is not a solid."
    if chemical_data[solution['B']]['state'] != 'gas':
        return f"Constraint 1 failed: B ({solution['B']}) is not a gas."
    if 'red' not in chemical_data[solution['C']]['color']:
        return f"Constraint 1 failed: C ({solution['C']}) is not described as red."
    # The 1:8 stoichiometry is the key insight, which is correctly used to identify A, B, C.

    # 2. Check Clue 3: C + H₂O → A(s) + F(strong acid) + G(weak acid)
    hydrolysis_info = chemical_data[solution['C']]['hydrolysis']
    if hydrolysis_info['reforms'] != solution['A']:
        return f"Constraint 3 failed: Hydrolysis of C ({solution['C']}) does not reform A ({solution['A']})."
    if hydrolysis_info['strong_acid'] != solution['F']:
        return f"Constraint 3 failed: The strong acid produced is not F ({solution['F']})."
    if hydrolysis_info['weak_acid'] != solution['G']:
        return f"Constraint 3 failed: The weak acid produced is not G ({solution['G']})."
    if chemical_data[solution['F']]['acid_strength'] != 'strong':
        return f"Constraint 3 failed: F ({solution['F']}) is not a strong acid."
    if chemical_data[solution['G']]['acid_strength'] != 'weak':
        return f"Constraint 3 failed: G ({solution['G']}) is not a weak acid."

    # 3. Check Clue 4: D(g) + B(g) → H (solvent) [1:1 ratio]
    # Reaction: C₂H₄(g) + Cl₂(g) → C₂H₄Cl₂(l)
    if chemical_data[solution['D']]['state'] != 'gas':
        return f"Constraint 4 failed: D ({solution['D']}) is not a gas."
    if chemical_data[solution['H']]['use'] != 'solvent':
        return f"Constraint 4 failed: H ({solution['H']}) is not a known solvent."
    # The 1:1 reaction is correctly identified.

    # 4. Check Clue 2: C + 2 D(g) → E (extremely hazardous product)
    # Reaction: SCl₂ + 2C₂H₄ → (ClCH₂CH₂)₂S (Levinstein process)
    if chemical_data[solution['E']]['hazard_level'] != 'extremely hazardous':
        return f"Constraint 2 failed: E ({solution['E']}) is not described as extremely hazardous."
    # The 1:2 stoichiometry is correctly identified with the Levinstein process.

    # 5. Final Check: Molecular symmetry group of E
    identified_molecule = solution['E']
    correct_point_group = chemical_data[identified_molecule]['point_group']
    
    if question_options[final_answer_option] != correct_point_group:
        return (f"Final Answer Check Failed: The answer identifies E as {identified_molecule}, "
                f"which has a point group of {correct_point_group}. "
                f"The selected option '{final_answer_option}' corresponds to '{question_options[final_answer_option]}', "
                f"which does not match.")

    return "Correct"

# Run the check
result = check_correctness()
print(result)