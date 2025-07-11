def propose_methanol_charges():
    """
    This function proposes a reasonable set of partial charges for methanol,
    prints them, and verifies that the molecule is charge-neutral.
    """

    # Define the charges for each atom based on the most reasonable choice
    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,
        "Methyl H2": 0.0400,
        "Methyl H3": 0.0400,
        "Oxygen": -0.6830,
        "Hydroxyl hydrogen": 0.4180
    }

    # Extract individual charges for clarity in the calculation
    c_charge = charges["Carbon"]
    h_charge = charges["Methyl H1"]  # All methyl hydrogens are the same
    o_charge = charges["Oxygen"]
    ho_charge = charges["Hydroxyl hydrogen"]

    # Calculate the total charge
    total_charge = c_charge + h_charge + h_charge + h_charge + o_charge + ho_charge

    # Print the proposed charges
    print("Proposed partial charges for methanol (CH3OH) in fundamental charge units (e):")
    for atom, charge in charges.items():
        # The :<20 format specifier left-aligns the string in a 20-character space
        # The : 8.4f format specifier formats the number as a float with 4 decimal places in a field of 8 characters, with a space for the sign
        print(f"{atom:<20}: {charge: 8.4f}")

    # Verify charge neutrality and print the equation
    print("\nVerification of charge neutrality:")
    print("The sum of charges for a neutral molecule must be zero.")
    
    # Print the full equation showing each number
    print(f"Calculation: {c_charge} + {h_charge} + {h_charge} + {h_charge} + {o_charge} + {ho_charge} = {total_charge:.4f}")
    
    if abs(total_charge) < 1e-9:
        print("The total charge is zero. The charge set is valid.")
    else:
        print(f"Warning: The total charge is {total_charge:.4f}, not zero. The charge set is invalid.")

# Execute the function
propose_methanol_charges()
