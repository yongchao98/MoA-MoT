import math

def evaluate_charge_sets():
    """
    Evaluates different sets of partial charges for methanol based on
    charge neutrality and chemical reasonableness.
    """
    choices = {
        'A': {'C': 0.5850, 'H_methyl': [-0.0860, -0.0860, -0.0860], 'O': -0.7597, 'H_hydroxyl': 0.4327},
        'B': {'C': 0.5845, 'H_methyl': [-0.0850, -0.0850, -0.0850], 'O': -0.7606, 'H_hydroxyl': 0.4332},
        'C': {'C': 0.5845, 'H_methyl': [-0.0857, -0.0857, -0.0857], 'O': -0.7606, 'H_hydroxyl': 0.4300},
        'D': {'C': 0.1450, 'H_methyl': [0.0400, 0.0400, 0.0400], 'O': -0.6830, 'H_hydroxyl': 0.4180},
        'E': {'C': 0.1450, 'H_methyl': [0.0400, 0.0300, 0.0300], 'O': -0.6830, 'H_hydroxyl': 0.4380}
    }

    best_choice = None
    best_score = -1

    print("Analyzing proposed partial charge sets for methanol (CH3OH):\n")

    for key, charges in choices.items():
        total_charge = charges['C'] + sum(charges['H_methyl']) + charges['O'] + charges['H_hydroxyl']
        
        # Criteria for a good model
        is_neutral = math.isclose(total_charge, 0.0, abs_tol=1e-4)
        methyl_h_symmetric = len(set(charges['H_methyl'])) == 1
        methyl_h_positive = all(h > 0 for h in charges['H_methyl'])
        
        print(f"--- Evaluating Choice {key} ---")
        print(f"Total charge: {total_charge:.4f}. Neutral: {is_neutral}")
        
        score = 0
        reasons = []
        if not is_neutral:
            reasons.append("Fails charge neutrality.")
        else:
            score += 1 # Base score for being neutral
            if not methyl_h_symmetric:
                reasons.append("Breaks methyl group symmetry, which is unreasonable for a standard model.")
            else:
                score += 1 # Bonus for symmetry
            
            if not methyl_h_positive:
                reasons.append("Assigns negative charges to methyl hydrogens, which is chemically questionable.")
            else:
                score += 1 # Bonus for chemical sense
        
        if not reasons:
            reasons.append("This is a reasonable charge set.")
            
        print(f"Analysis: {' '.join(reasons)}\n")

        if score > best_score:
            best_score = score
            best_choice = key

    print("--- Conclusion ---")
    print(f"Choice {best_choice} is the most reasonable proposal.")
    print("It is charge-neutral, maintains the symmetry of the methyl hydrogens, and assigns charges consistent with electronegativity principles (positive methyl hydrogens).\n")
    
    # Print the final proposed charges and the neutrality equation
    final_charges = choices[best_choice]
    c_charge = final_charges['C']
    h_methyl_charges = final_charges['H_methyl']
    o_charge = final_charges['O']
    h_hydroxyl_charge = final_charges['H_hydroxyl']
    total = c_charge + sum(h_methyl_charges) + o_charge + h_hydroxyl_charge

    print("Proposed Partial Charges:")
    print(f"Carbon: \t{c_charge:.4f}")
    print(f"Methyl H1:\t{h_methyl_charges[0]:.4f}")
    print(f"Methyl H2:\t{h_methyl_charges[1]:.4f}")
    print(f"Methyl H3:\t{h_methyl_charges[2]:.4f}")
    print(f"Oxygen: \t{o_charge:.4f}")
    print(f"Hydroxyl H:\t{h_hydroxyl_charge:.4f}")
    print("\nCharge Neutrality Check:")
    print(f"Carbon ({c_charge:+.4f}) + Methyl H1 ({h_methyl_charges[0]:+.4f}) + Methyl H2 ({h_methyl_charges[1]:+.4f}) + Methyl H3 ({h_methyl_charges[2]:+.4f}) + Oxygen ({o_charge:+.4f}) + Hydroxyl H ({h_hydroxyl_charge:+.4f}) = {total:.4f}")

evaluate_charge_sets()
print("\n<<<D>>>")