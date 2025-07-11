import sys

def propose_methanol_charges():
    """
    Proposes and prints a reasonable set of partial charges for an all-atom
    model of methanol (CH3OH) and verifies its neutrality.
    """
    # Based on analysis, this set is chemically reasonable, symmetric, and neutral.
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl hydrogen': 0.4180
    }

    print("Proposed partial charge assignments for Methanol (CH3OH):")
    print("Atom\t\t\tCharge (e)")
    print("---------------------------------------")
    for atom, charge in charges.items():
        # Using a fixed width for alignment
        print(f"{atom:<20}\t{charge:+.4f}")

    # Verify that the total charge is neutral and show the equation
    total_charge = sum(charges.values())
    charge_values = list(charges.values())
    
    # Building the string for the summation equation.
    # We will use ' ' for positive numbers to align with '-' for negative numbers.
    equation_parts = []
    for charge in charge_values:
        if charge >= 0:
            equation_parts.append(f" {charge:.4f}")
        else:
            equation_parts.append(f"{charge:.4f}")

    # Joining with ' +' to create the sum expression.
    equation_str = " +".join(equation_parts)
    
    print("\nVerification of Net Charge Neutrality:")
    print(f"Sum = {equation_str} = {total_charge:+.4f}")

    # A small tolerance is used for floating point arithmetic comparison
    if abs(total_charge) < 1e-9:
        print("Conclusion: The total charge is 0. The molecule is electrically neutral.")
    else:
        # This part will not be reached for the chosen values but is good practice.
        print(f"Conclusion: Warning! The total charge is {total_charge:f}. The molecule is not neutral.")


if __name__ == '__main__':
    propose_methanol_charges()