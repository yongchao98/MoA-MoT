import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for methanol (CH3OH)
    and verifies that the molecule is electrically neutral.
    """

    # Partial charges based on the most chemically sound option.
    charge_C = 0.1450
    charge_H_methyl = 0.0400
    charge_O = -0.6830
    charge_H_hydroxyl = 0.4180

    print("Proposed Partial Charges for Methanol (CH3OH):")
    print("-" * 45)
    # The f-string formatting ensures trailing zeros are kept for precision
    print(f"Carbon (C):           {charge_C: 9.4f}")
    print(f"Methyl Hydrogen (H):  {charge_H_methyl: 9.4f} (for each of 3)")
    print(f"Oxygen (O):           {charge_O: 9.4f}")
    print(f"Hydroxyl Hydrogen (H):{charge_H_hydroxyl: 9.4f}")
    print("-" * 45)

    # Verify that the sum of charges is zero.
    total_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl

    # Use a small tolerance for floating point comparison
    if abs(total_charge) > 1e-9:
        print(f"\nWarning: The total charge is {total_charge:.4f}, which is not neutral.", file=sys.stderr)
    else:
        # Round the result to avoid floating point artifacts like -0.0
        total_charge = 0.0

    print("\nVerifying charge neutrality of the molecule:")
    # Print the equation with all the numbers, as requested.
    print(f"({charge_C:.4f}) + 3*({charge_H_methyl:.4f}) + ({charge_O:.4f}) + ({charge_H_hydroxyl:.4f}) = {total_charge:.4f}")
    print("The molecule is electrically neutral.")

if __name__ == "__main__":
    propose_methanol_charges()
