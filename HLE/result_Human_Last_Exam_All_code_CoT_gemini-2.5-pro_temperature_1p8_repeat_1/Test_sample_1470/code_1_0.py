import sys

def evaluate_methanol_charges():
    """
    This function evaluates proposed partial charge assignments for methanol.
    It checks for charge neutrality and chemical reasonableness for each choice.
    """

    choices = {
        'A': {
            'C': 0.5850, 'H1': -0.0860, 'H2': -0.0860, 'H3': -0.0860,
            'O': -0.7597, 'H_ox': 0.4327,
            'symmetric_H': True
        },
        'B': {
            'C': 0.5845, 'H1': -0.0850, 'H2': -0.0850, 'H3': -0.0850,
            'O': -0.7606, 'H_ox': 0.4332,
            'symmetric_H': True
        },
        'C': {
            'C': 0.5845, 'H1': -0.0857, 'H2': -0.0857, 'H3': -0.0857,
            'O': -0.7606, 'H_ox': 0.4300,
            'symmetric_H': True
        },
        'D': {
            'C': 0.1450, 'H1': 0.0400, 'H2': 0.0400, 'H3': 0.0400,
            'O': -0.6830, 'H_ox': 0.4180,
            'symmetric_H': True
        },
        'E': {
            'C': 0.1450, 'H1': 0.0400, 'H2': 0.0300, 'H3': 0.0300,
            'O': -0.6830, 'H_ox': 0.4380,
            'symmetric_H': False
        }
    }

    print("Evaluating Proposed Partial Charges for Methanol (CH3OH)\n")
    best_choice = None
    
    for choice, charges in choices.items():
        c, h1, h2, h3, o, h_ox = charges['C'], charges['H1'], charges['H2'], charges['H3'], charges['O'], charges['H_ox']
        total_charge = c + h1 + h2 + h3 + o + h_ox
        
        print(f"--- Analyzing Choice {choice} ---")
        print(f"Proposed Charges:")
        print(f"Carbon\t{c:.4f}")
        print(f"Methyl H1:\t{h1:.4f}")
        print(f"Methyl H2:\t{h2:.4f}")
        print(f"Methyl H3:\t{h3:.4f}")
        print(f"Oxygen:\t{o:.4f}")
        print(f"Hydroxyl H:\t{h_ox:.4f}\n")
        
        print("Summing charges to check for neutrality:")
        print(f"Equation: {c:.4f} (C) + {h1:.4f} (H1) + {h2:.4f} (H2) + {h3:.4f} (H3) + {o:.4f} (O) + {h_ox:.4f} (H_ox)")
        # Use a high precision for the sum to avoid floating point artifacts
        print(f"Total Charge = {total_charge:.4f}\n")

        # Analysis
        is_neutral = abs(total_charge) < 1e-9
        is_chemically_sound = h1 > 0 and h2 > 0 and h3 > 0 and o < 0 and h_ox > 0
        is_symmetric = charges['symmetric_H']

        print("Evaluation:")
        if not is_neutral:
            print("Verdict: INVALID. The molecule is not charge-neutral.")
        elif not is_chemically_sound:
            print("Verdict: INVALID. Charges are not chemically reasonable (e.g., methyl hydrogens should be positive).")
        elif not is_symmetric:
             print("Verdict: NOT IDEAL. Breaks the symmetry of the methyl group, which is unusual for a general model.")
        else:
            print("Verdict: VALID. The molecule is charge-neutral and charges are chemically reasonable.")
            if best_choice is None: # Found the first valid choice
                best_choice = choice
        print("-" * 25 + "\n")

    print("\n--- Final Conclusion ---")
    print("Based on the analysis, only one option meets all the criteria for a reasonable molecular model:")
    print(f"Choice {best_choice} is the only one that is both charge-neutral and has chemically plausible charge assignments for all atoms, while respecting the molecular symmetry.")
    print("\n- Choices A, B, and C propose negative charges for the methyl hydrogens, which contradicts chemical intuition based on electronegativity.")
    print("- Choices B and C also fail the charge neutrality test.")
    print("- Choice E is charge-neutral, but it breaks the chemical equivalence of the three methyl hydrogens, which is not a good starting point for a general model.")
    print(f"\nTherefore, Choice {best_choice} provides the most reasonable set of partial charges.")


if __name__ == '__main__':
    evaluate_methanol_charges()