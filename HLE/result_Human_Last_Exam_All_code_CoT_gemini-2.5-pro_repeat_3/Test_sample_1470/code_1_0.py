def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for a methanol
    molecule and verifies that the molecule is charge-neutral.
    """
    
    # Proposed partial charges for Methanol (CH3OH)
    # Based on chemical principles and common force field models.
    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,
        "Methyl H2": 0.0400,
        "Methyl H3": 0.0400,
        "Oxygen": -0.6830,
        "Hydroxyl hydrogen": 0.4180
    }
    
    # Calculate the total charge
    total_charge = sum(charges.values())
    
    # Print the proposed charges
    print("Proposed Partial Charges for Methanol (in elementary charge units):\n")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")
        
    # Verify and print the total charge
    print("\n------------------------------------")
    print(f"Verification of Charge Neutrality:")
    # The final equation is the sum of all charges.
    charge_sum_str = " + ".join([f"{v:.4f}" for v in charges.values()])
    print(f"Sum of charges = {charge_sum_str} = {total_charge:.4f}")
    
    if abs(total_charge) < 1e-9:
        print("The molecule is charge-neutral. These are reasonable charge assignments.")
    else:
        print(f"Warning: The molecule is not charge-neutral (Total charge: {total_charge:.4f}).")

propose_methanol_charges()
<<<D>>>