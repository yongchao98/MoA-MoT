import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for methanol
    by evaluating the given options against chemical and physical principles.

    The chosen set of charges must:
    1. Sum to zero for a neutral molecule.
    2. Respect the chemical symmetry of the methyl group (all three H's have the same charge).
    3. Be chemically reasonable based on electronegativity.
    """

    # Charges from option D, which is the most reasonable choice.
    # O is most electronegative -> most negative charge.
    # H bonded to O is highly polarized -> most positive charge.
    # C is bonded to O -> positive charge.
    # Methyl H's are slightly positive.
    # The three methyl H's are symmetric.
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl H': 0.4180
    }

    # Verify that the total charge is zero.
    c_charge = charges['Carbon']
    h_methyl_charge = charges['Methyl H1'] # They are all the same
    o_charge = charges['Oxygen']
    h_hydroxyl_charge = charges['Hydroxyl H']

    total_charge = c_charge + 3 * h_methyl_charge + o_charge + h_hydroxyl_charge

    # Print the proposed charges
    print("Proposed Partial Charges for Methanol (CH3OH):")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")

    print("\nVerification of Net Charge Neutrality:")
    # sys.stdout.write is used to avoid adding a space between elements that print() does
    sys.stdout.write(f"{c_charge:.4f} (C) + 3 * {h_methyl_charge:.4f} (H_methyl) + ({o_charge:.4f}) (O) + {h_hydroxyl_charge:.4f} (H_hydroxyl) = {total_charge:.4f}\n")
    
    if abs(total_charge) < 1e-9:
        print("\nThe sum of the partial charges is zero. This is a valid charge set for a neutral molecule.")
    else:
        print(f"\nWarning: The sum of the partial charges is {total_charge:.4f}, not zero.")

propose_methanol_charges()