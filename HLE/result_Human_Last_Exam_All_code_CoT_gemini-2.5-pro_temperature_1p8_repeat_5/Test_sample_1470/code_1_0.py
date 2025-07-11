def assign_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for methanol,
    verifies the net charge is zero, and prints the assignments.
    """
    # Proposed partial charges from the most reasonable choice
    charge_C = 0.1450
    charge_H_methyl = 0.0400
    charge_O = -0.6830
    charge_H_hydroxyl = 0.4180
    
    # Store atom charges in a dictionary for clear labeling
    charges = {
        'Carbon': charge_C,
        'Methyl H1': charge_H_methyl,
        'Methyl H2': charge_H_methyl,
        'Methyl H3': charge_H_methyl,
        'Oxygen': charge_O,
        'Hydroxyl Hydrogen': charge_H_hydroxyl
    }
    
    print("Proposed Partial Charges for Methanol (CH3OH):")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")
        
    # Verify that the total charge sums to zero
    total_charge = charge_C + 3 * charge_H_methyl + charge_O + charge_H_hydroxyl
    
    print("\nVerification of Net Charge Neutrality:")
    print(f"Equation: C + 3*(H_methyl) + O + H_hydroxyl = Total Charge")
    # The prompt requests each number in the final equation
    print(f"Calculation: {charge_C:.4f} + 3*({charge_H_methyl:.4f}) + ({charge_O:.4f}) + {charge_H_hydroxyl:.4f} = {total_charge:.4f}")
    if round(total_charge, 4) == 0.0:
        print("The molecule is electrically neutral.")
    else:
        print(f"Warning: The molecule is not neutral. Total charge is {total_charge:.4f}")

assign_methanol_charges()
<<<D>>>