def check_charges():
    """
    Proposes and verifies a reasonable set of partial charges for methanol.
    
    This function selects the most chemically sound and charge-neutral option
    from a predefined list and prints the charge assignment along with a
    verification of its neutrality.
    """
    
    # Methanol atoms: 1 Carbon, 3 methyl Hydrogens, 1 Oxygen, 1 hydroxyl Hydrogen
    
    # Based on chemical principles and charge neutrality, Option D is the best choice.
    # Charges from Option D
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl H': 0.4180
    }
    
    print("Proposed Partial Charges for Methanol (CH3OH) Model:")
    print("--------------------------------------------------")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f} e")
        
    # Verify charge neutrality
    total_charge = sum(charges.values())
    
    print("\nVerifying charge neutrality of the proposed model:")
    
    # Building the equation string with all numbers
    charge_values = [v for v in charges.values()]
    equation_str = " + ".join([f"({c})" if c < 0 else str(c) for c in charge_values])
    
    print(f"Sum = {equation_str} = {total_charge:.4f}")
    
    if abs(total_charge) < 1e-9:
        print("\nThe total charge is 0. This is a valid and chemically reasonable set of partial charges.")
    else:
        print(f"\nWarning: The total charge is {total_charge:.4f}, not 0. This set of charges is not valid.")

if __name__ == "__main__":
    check_charges()