def propose_methanol_charges():
    """
    Proposes a set of reasonable partial charges for a methanol (CH3OH)
    molecule and verifies that the molecule is electroneutral.
    """
    # These charges are based on established principles of electronegativity and symmetry.
    # Among the given options, this set is the only one that is neutral,
    # symmetric for the methyl group, and chemically sound.
    # C: Carbon, Hc: Methyl Hydrogen, O: Oxygen, Ho: Hydroxyl Hydrogen
    charge_C = 0.1450
    charge_Hc = 0.0400
    charge_O = -0.6830
    charge_Ho = 0.4180

    print("Proposed reasonable partial charges for methanol (CH3OH) in fundamental charge units:")
    print("-" * 40)
    print(f"Carbon:\t\t\t{charge_C:.4f}")
    # The three methyl hydrogens are chemically equivalent
    print(f"Methyl H1:\t\t{charge_Hc:.4f}")
    print(f"Methyl H2:\t\t{charge_Hc:.4f}")
    print(f"Methyl H3:\t\t{charge_Hc:.4f}")
    print(f"Oxygen:\t\t\t{charge_O:.4f}")
    print(f"Hydroxyl hydrogen:\t{charge_Ho:.4f}")
    print("-" * 40)
    
    # Verify the electroneutrality of the molecule by summing the charges
    total_charge = charge_C + 3 * charge_Hc + charge_O + charge_Ho

    print("Verification of electroneutrality:")
    print("Total Charge = (C) + 3*(Methyl H) + (O) + (Hydroxyl H)")
    
    # Print the equation with the actual numbers, as requested.
    print(f"Total Charge = {charge_C:.4f} + 3*({charge_Hc:.4f}) + ({charge_O:.4f}) + {charge_Ho:.4f}")
    
    # To be even more explicit, we can show the result of the multiplication.
    print(f"Total Charge = {charge_C:.4f} + {3 * charge_Hc:.4f} + {charge_O:.4f} + {charge_Ho:.4f}")
    
    # Print the final sum
    print(f"Final net charge: {total_charge:.4f}")

propose_methanol_charges()