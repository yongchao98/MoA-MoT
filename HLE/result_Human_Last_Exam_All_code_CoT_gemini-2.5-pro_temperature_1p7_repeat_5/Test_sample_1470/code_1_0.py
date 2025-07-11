def propose_methanol_charges():
    """
    Proposes a reasonable set of partial charges for methanol and verifies neutrality.
    """

    # The chosen set of charges is from Choice D, as it is the most chemically
    # and physically reasonable option provided.
    charges = {
        "Carbon": 0.1450,
        "Methyl H1": 0.0400,
        "Methyl H2": 0.0400,
        "Methyl H3": 0.0400,
        "Oxygen": -0.6830,
        "Hydroxyl hydrogen": 0.4180
    }

    # Print the proposed charges for the user.
    print("Proposed reasonable partial charges for methanol (in fundamental charge units):\n")
    for atom, charge in charges.items():
        print(f"{atom:<20}: {charge:>9.4f}")

    # Verify that the total charge of the molecule is zero.
    print("\nVerification that the net charge is zero:")
    
    charge_values = list(charges.values())
    total_charge = sum(charge_values)

    # Format each number to always show its sign and 4 decimal places for the equation.
    equation_parts = [f"({charge:+.4f})" for charge in charge_values]
    equation_string = " + ".join(equation_parts)

    print(f"{equation_string} = {total_charge:.4f}")

propose_methanol_charges()
<<<D>>>