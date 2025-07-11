import math

def propose_and_verify_charges():
    """
    Proposes a set of reasonable partial charges for methanol and verifies
    that the molecule is charge-neutral.
    """
    # Proposing a charge set inspired by successful existing force fields.
    # This set is chemically reasonable, symmetric for the methyl group,
    # and maintains the overall neutrality of the molecule.
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl hydrogen': 0.4180
    }

    print("Proposed partial charge assignments for Methanol (CH3OH):")
    for atom, charge in charges.items():
        # Ensure consistent formatting
        print(f"{atom}:\t{charge:.4f}")

    # Verify charge neutrality
    total_charge = sum(charges.values())

    print("\nVerifying the total charge of the molecule:")
    
    # Building the equation string with each number
    charge_values = [
        charges['Carbon'], 
        charges['Methyl H1'], 
        charges['Methyl H2'], 
        charges['Methyl H3'], 
        charges['Oxygen'], 
        charges['Hydroxyl hydrogen']
    ]
    
    # Join the formatted numbers with ' + '
    # For negative numbers, this will look like "... + -0.6830 + ..." which is arithmetically clear
    equation = " + ".join(f"{c:.4f}" for c in charge_values)
    
    print(f"{equation} = {total_charge:.4f}")

    # Check if the sum is effectively zero
    if math.isclose(total_charge, 0.0, abs_tol=1e-5):
        print("\nThe sum of the partial charges is zero.")
        print("This is a valid and reasonable set of charges for a simulation.")
    else:
        print("\nWarning: The sum of the partial charges is not zero.")
        print("This set of charges is not valid for a neutral molecule.")

if __name__ == "__main__":
    propose_and_verify_charges()