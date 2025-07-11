import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for methanol,
    prints them, and verifies that the molecule is charge-neutral.
    """
    
    # Based on analysis, choice D is the only one that is charge-neutral,
    # respects chemical symmetry, and follows electronegativity principles.
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl hydrogen': 0.4180
    }

    # Extracting charge values for the calculation
    q_C = charges['Carbon']
    # All methyl hydrogens are the same
    q_Hc = charges['Methyl H1'] 
    q_O = charges['Oxygen']
    q_Ho = charges['Hydroxyl hydrogen']

    # Calculate the total charge of the molecule
    total_charge = q_C + 3 * q_Hc + q_O + q_Ho

    print("Proposed reasonable partial charge assignments for methanol (in fundamental charge units):")
    # Using a format specifier to ensure consistent decimal places
    for atom, charge in charges.items():
        print(f"{atom:<20}{charge: .4f}")

    print("\nTo be a valid model, the sum of partial charges must be zero.")
    print("Verification of charge neutrality:")
    # The final print statement explicitly shows each number in the equation.
    print(f"Sum = {q_C:.4f} (C) + 3 * {q_Hc:.4f} (H_methyl) + {q_O:.4f} (O) + {q_Ho:.4f} (H_hydroxyl) = {total_charge:.4f}")

if __name__ == '__main__':
    propose_methanol_charges()
