import math

def propose_methanol_charges():
    """
    This function defines and verifies a reasonable set of partial charges for methanol.
    It selects the only chemically and mathematically sound option from the choices provided.
    """
    
    # Partial charges from the most reasonable choice (D)
    charge_C = 0.1450
    charge_H_methyl = 0.0400  # All methyl hydrogens are equivalent
    charge_O = -0.6830
    charge_H_hydroxyl = 0.4180
    
    num_methyl_H = 3
    
    # Calculate the total charge to verify neutrality
    total_charge = charge_C + (num_methyl_H * charge_H_methyl) + charge_O + charge_H_hydroxyl
    
    print("Proposed Reasonable Partial Charges for Methanol (CH3OH):")
    print(f"{'Atom':<20} {'Charge (e)':<10}")
    print("-" * 30)
    print(f"{'Carbon':<20} {charge_C:>10.4f}")
    print(f"{'Methyl H1':<20} {charge_H_methyl:>10.4f}")
    print(f"{'Methyl H2':<20} {charge_H_methyl:>10.4f}")
    print(f"{'Methyl H3':<20} {charge_H_methyl:>10.4f}")
    print(f"{'Oxygen':<20} {charge_O:>10.4f}")
    print(f"{'Hydroxyl hydrogen':<20} {charge_H_hydroxyl:>10.4f}")
    print("-" * 30)
    
    print("\nVerification of Net Charge:")
    # Using 'math.copysign' to correctly display the signs for the equation
    print(f"{charge_C:.4f} (C) + {num_methyl_H}*({charge_H_methyl:.4f}) (methyl H's) + ({charge_O:.4f}) (O) + {charge_H_hydroxyl:.4f} (hydroxyl H) = {total_charge:.4f}")

if __name__ == '__main__':
    propose_methanol_charges()