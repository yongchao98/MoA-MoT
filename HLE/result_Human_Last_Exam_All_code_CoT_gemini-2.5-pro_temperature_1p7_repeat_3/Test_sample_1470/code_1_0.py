import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for a new all-atom
    model of methanol (CH3OH) and verifies its neutrality.
    The selected charges are based on physical and chemical principles,
    respecting symmetry and electronegativity.
    """
    
    # Based on evaluation, the charges from option D are the most reasonable.
    # Structure:
    #      H
    #      |
    #  H - C - O - H
    #      |
    #      H

    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,  # All three methyl hydrogens are equivalent
        "Oxygen": -0.6830,
        "Hydroxyl H": 0.4180
    }

    c_charge = charges["Carbon"]
    h_methyl_charge = charges["Methyl H1"]
    o_charge = charges["Oxygen"]
    h_hydroxyl_charge = charges["Hydroxyl H"]

    print("Proposed Partial Charges for Methanol (in fundamental charge units 'e'):")
    print("----------------------------------------------------------------------")
    # Using format specifiers for clean alignment and to show the sign.
    print(f"Atom               | Charge")
    print("------------------ | ---------")
    print(f"Carbon             | {c_charge:+.4f}")
    print(f"Methyl H1          | {h_methyl_charge:+.4f}")
    print(f"Methyl H2          | {h_methyl_charge:+.4f}")
    print(f"Methyl H3          | {h_methyl_charge:+.4f}")
    print(f"Oxygen             | {o_charge:+.4f}")
    print(f"Hydroxyl Hydrogen  | {h_hydroxyl_charge:+.4f}")
    print("----------------------------------------------------------------------")

    # Verify that the total charge of the molecule is zero.
    # This is a critical check for any charge model of a neutral molecule.
    total_charge = c_charge + 3 * h_methyl_charge + o_charge + h_hydroxyl_charge
    
    print("\nVerification of Net Molecular Charge:")
    # Print each number in the final equation as requested.
    equation_str = (
        f"({c_charge:+.4f}) + 3*({h_methyl_charge:+.4f}) + "
        f"({o_charge:+.4f}) + ({h_hydroxyl_charge:+.4f})"
    )
    print(f"Summation: {equation_str} = {total_charge:.4f}")
    
    if abs(total_charge) < 1e-9:
        print("Result: The total charge is zero, as required for a neutral molecule.")
    else:
        print(f"Warning: The total charge is {total_charge:.4f}, not zero.", file=sys.stderr)

if __name__ == '__main__':
    propose_methanol_charges()