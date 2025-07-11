import math

def evaluate_methanol_charges():
    """
    Evaluates different sets of partial charges for methanol based on
    net neutrality, symmetry, and chemical reasonableness.
    """
    # Define the charge sets from the answer choices
    # Format: [C, H_methyl1, H_methyl2, H_methyl3, O, H_hydroxyl]
    charge_sets = {
        'A': [0.5850, -0.0860, -0.0860, -0.0860, -0.7597, 0.4327],
        'B': [0.5845, -0.0850, -0.0850, -0.0850, -0.7606, 0.4332],
        'C': [0.5845, -0.0857, -0.0857, -0.0857, -0.7606, 0.4300],
        'D': [0.1450, 0.0400, 0.0400, 0.0400, -0.6830, 0.4180],
        'E': [0.1450, 0.0400, 0.0300, 0.0300, -0.6830, 0.4380]
    }

    print("--- Analysis of Partial Charge Sets for Methanol ---")
    best_option = None
    
    for option, charges in charge_sets.items():
        total_charge = sum(charges)
        # Check for symmetry in methyl hydrogens
        h_methyl_symmetric = (charges[1] == charges[2] == charges[3])
        # Check for chemically sound charges (O negative, H_methyl positive)
        chemically_sound = (charges[4] < 0 and charges[1] > 0)

        print(f"\nEvaluating Option {option}:")
        print(f"Charges: C={charges[0]}, H_m={charges[1]}/{charges[2]}/{charges[3]}, O={charges[4]}, H_o={charges[5]}")
        print(f"Total Charge: {total_charge:.4f}")

        # Check conditions
        if not math.isclose(total_charge, 0.0, abs_tol=1e-4):
            print("Result: Fails. The molecule is not neutral.")
        elif not h_methyl_symmetric:
            print("Result: Fails. Methyl hydrogens are not symmetrical, which is physically unreasonable.")
        elif not chemically_sound:
            print("Result: Plausible, but less common. Methyl hydrogens have non-positive charges.")
        else:
            print("Result: Pass. The set is neutral, symmetrical, and chemically sound.")
            best_option = option

    print("\n--- Conclusion ---")
    if best_option:
        print(f"Option {best_option} is the most reasonable choice.")
        print("It maintains net neutrality, treats equivalent hydrogens symmetrically, and aligns with chemical principles of electronegativity.")
        
        chosen_charges = charge_sets[best_option]
        c, h1, h2, h3, o, ho = chosen_charges
        
        print("\nProposed Partial Charge Assignment (Option D):")
        print(f"Carbon: {c:.4f}")
        print(f"Methyl H1: {h1:.4f}")
        print(f"Methyl H2: {h2:.4f}")
        print(f"Methyl H3: {h3:.4f}")
        print(f"Oxygen: {o:.4f}")
        print(f"Hydroxyl hydrogen: {ho:.4f}")
        
        print("\nVerification of Net Charge Equation:")
        # Output each number in the final equation
        print(f"{c:.4f} (C) + {h1:.4f} (H1) + {h2:.4f} (H2) + {h3:.4f} (H3) + {o:.4f} (O) + {ho:.4f} (Ho) = {sum(chosen_charges):.4f}")
    else:
        print("None of the options perfectly satisfy all criteria for a standard model.")

evaluate_methanol_charges()