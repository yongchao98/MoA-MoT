import collections

def propose_methanol_charges():
    """
    Proposes a reasonable set of partial charges for methanol,
    prints each charge, and verifies that the molecule is charge-neutral.
    """
    # A reasonable set of partial charges for the atoms in methanol (CH3OH)
    # C, H(methyl), H(methyl), H(methyl), O, H(hydroxyl)
    charges = collections.OrderedDict([
        ("Carbon", 0.1450),
        ("Methyl H1", 0.0400),
        ("Methyl H2", 0.0400),
        ("Methyl H3", 0.0400),
        ("Oxygen", -0.6830),
        ("Hydroxyl hydrogen", 0.4180)
    ])

    print("Proposed Partial Charge Assignments for Methanol:")
    print("---------------------------------------------")
    
    charge_values = list(charges.values())
    total_charge = sum(charge_values)

    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")
    
    print("\n---------------------------------------------")
    print("Verification of Charge Neutrality:")
    
    # Building and printing the full equation for the sum
    equation_parts = []
    for charge in charge_values:
        if charge < 0:
            equation_parts.append(f"({charge:.4f})")
        else:
            equation_parts.append(f"{charge:.4f}")
    
    equation_str = " + ".join(equation_parts)
    print(f"Sum of charges = {equation_str} = {total_charge:.4f}")

    if abs(total_charge) < 1e-6:
        print("\nThe proposed charges result in a neutral molecule.")
    else:
        print(f"\nWarning: The molecule is not charge neutral. Total charge is {total_charge:.4f}")

propose_methanol_charges()