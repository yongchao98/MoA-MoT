import math

def analyze_methanol_charges():
    """
    Analyzes potential partial charge models for methanol (CH3OH)
    and recommends the most reasonable set based on fundamental principles.
    """
    # Structure of methanol: 1 Carbon, 3 Methyl Hydrogens, 1 Oxygen, 1 Hydroxyl Hydrogen
    choices = {
        'A': {'C': 0.5850, 'HC': -0.0860, 'O': -0.7597, 'HO': 0.4327},
        'B': {'C': 0.5845, 'HC': -0.0850, 'O': -0.7606, 'HO': 0.4332},
        'C': {'C': 0.5845, 'HC_list': [-0.0857, -0.0857, -0.0857], 'O': -0.7606, 'HO': 0.4300},
        'D': {'C': 0.1450, 'HC': 0.0400, 'O': -0.6830, 'HO': 0.4180},
        'E': {'C': 0.1450, 'HC_list': [0.0400, 0.0300, 0.0300], 'O': -0.6830, 'HO': 0.4380}
    }

    best_choice = None
    
    print("Evaluating Proposed Methanol Charge Models")
    print("="*45)

    for choice_id, charges in choices.items():
        is_symmetric = 'HC_list' not in charges
        
        # Unify the data structure for easier calculation
        if is_symmetric:
            c = charges['C']
            hc1 = hc2 = hc3 = charges['HC']
        else: # Handles asymmetric choices C and E
            c = charges['C']
            hc1, hc2, hc3 = charges['HC_list']

        o = charges['O']
        ho = charges['HO']
        
        total_charge = c + hc1 + hc2 + hc3 + o + ho
        
        # Check against criteria
        is_neutral = math.isclose(total_charge, 0, abs_tol=1e-4)
        is_chemically_sound = (hc1 > 0 and hc2 > 0 and hc3 > 0)

        # Verdict
        if is_neutral and is_symmetric and is_chemically_sound:
            best_choice = choice_id
            
    if best_choice:
        print(f"\nConclusion: Choice '{best_choice}' provides the most reasonable charge assignment.\n")
        print("It satisfies the key requirements for a molecular model:")
        print("1. Net charge is zero (neutrality).")
        print("2. Equivalent atoms have identical charges (symmetry).")
        print("3. Charges reflect chemical electronegativity (chemical soundness).\n")
        
        final_model = choices[best_choice]
        c, hc, o, ho = final_model['C'], final_model['HC'], final_model['O'], final_model['HO']
        total = c + (3 * hc) + o + ho
        
        print("Proposed Partial Charge Assignment:")
        print(f"{'Atom':<20} {'Charge (e)':<15}")
        print("-" * 35)
        print(f"{'Carbon':<20} {c:<+15.4f}")
        print(f"{'Methyl Hydrogen (x3)':<20} {hc:<+15.4f}")
        print(f"{'Oxygen':<20} {o:<+15.4f}")
        print(f"{'Hydroxyl Hydrogen':<20} {ho:<+15.4f}")
        print("-" * 35)

        print("\nVerification of Net Charge Neutrality:")
        print("Equation: Charge(C) + 3*Charge(H_methyl) + Charge(O) + Charge(H_hydroxyl) = Total")
        print(f"Values:   {c:.4f} + 3*({hc:.4f}) + ({o:.4f}) + {ho:.4f} = {total:.4f}")

    else:
        print("\nNo single choice perfectly satisfies all criteria for an ideal model.")


if __name__ == '__main__':
    analyze_methanol_charges()
