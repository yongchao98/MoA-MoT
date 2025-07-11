import collections

def propose_methanol_charges():
    """
    Proposes a reasonable set of partial charges for a methanol molecule,
    prints the assignments, and verifies charge neutrality.
    """
    # Using an ordered dictionary to maintain a logical atom order for printing.
    # These charges are from choice D, which is chemically reasonable and neutral.
    charges = collections.OrderedDict([
        ('Carbon', 0.1450),
        ('Methyl H1', 0.0400),
        ('Methyl H2', 0.0400),
        ('Methyl H3', 0.0400),
        ('Oxygen', -0.6830),
        ('Hydroxyl H', 0.4180)
    ])

    print("Proposed Partial Charges for Methanol (in fundamental charge units):\n")
    for atom, charge in charges.items():
        # Print each atom and its assigned charge, formatted for clarity
        print(f"{atom:<12}: {charge:8.4f}")

    # Verify that the sum of the charges is zero
    total_charge = sum(charges.values())

    # Build the string for the summation equation
    charge_values = list(charges.values())
    equation_str = " + ".join(f"({v:.4f})" for v in charge_values)

    print("\nVerification of charge neutrality:")
    # Using round() to handle potential floating point inaccuracies for the final print
    print(f"Sum = {equation_str} = {round(total_charge, 4)}")

propose_methanol_charges()