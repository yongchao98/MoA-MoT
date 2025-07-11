import sys

# This script proposes a reasonable set of partial charges for methanol (CH3OH)
# based on an evaluation of several options. The chosen set adheres to principles
# of charge neutrality, symmetry, and chemical plausibility.

def propose_methanol_charges():
    """
    Prints the chosen partial charges for each atom in methanol and verifies
    that the net charge of the molecule is zero.
    """
    
    # Based on the analysis, option D is the most reasonable set of charges.
    # C=0.1450, H(methyl)=0.0400, O=-0.6830, H(hydroxyl)=0.4180
    
    c_charge = 0.1450
    h_methyl_charge = 0.0400
    o_charge = -0.6830
    h_hydroxyl_charge = 0.4180
    
    print("Proposed partial charge assignments for a new all-atom model of methanol:")
    print("-" * 70)
    # Using python's string formatting to align the output nicely.
    print(f"{'Atom':<20} {'Partial Charge (e)':<20}")
    print("-" * 70)
    print(f"{'Carbon':<20} {c_charge:>20.4f}")
    print(f"{'Methyl H1':<20} {h_methyl_charge:>20.4f}")
    print(f"{'Methyl H2':<20} {h_methyl_charge:>20.4f}")
    print(f"{'Methyl H3':<20} {h_methyl_charge:>20.4f}")
    print(f"{'Oxygen':<20} {o_charge:>20.4f}")
    print(f"{'Hydroxyl hydrogen':<20} {h_hydroxyl_charge:>20.4f}")
    print("-" * 70)
    
    # Verify that the total charge sums to zero
    total_charge = c_charge + (3 * h_methyl_charge) + o_charge + h_hydroxyl_charge
    
    # Print the verification equation with all the numbers
    print("\nVerification of Net Charge:")
    print(f"Carbon({c_charge:.4f}) + 3*H_methyl({h_methyl_charge:.4f}) + Oxygen({o_charge:.4f}) + H_hydroxyl({h_hydroxyl_charge:.4f}) = {total_charge:.4f}")

if __name__ == '__main__':
    propose_methanol_charges()
