import collections

def propose_methanol_charges():
    """
    Proposes a set of reasonable partial charges for a new all-atom model
    of methanol and verifies the neutrality of the molecule.
    """

    # Proposed partial charges for methanol (CH3OH) based on chemical principles.
    # The three methyl hydrogens are treated as equivalent.
    charges = collections.OrderedDict([
        ("Carbon", 0.1450),
        ("Methyl H1", 0.0400),
        ("Methyl H2", 0.0400),
        ("Methyl H3", 0.0400),
        ("Oxygen", -0.6830),
        ("Hydroxyl hydrogen", 0.4180)
    ])

    print("Proposed partial charges for methanol (in elementary charge units):")
    print("-" * 60)

    # Print each atom's charge
    for atom, charge in charges.items():
        print(f"{atom:<20}: {charge:>9.4f}")

    print("-" * 60)

    # Verify that the total charge is zero
    total_charge = sum(charges.values())
    
    # Create the equation string
    charge_values = [f"({v:.4f})" if v < 0 else f"{v:.4f}" for v in charges.values()]
    equation_str = " + ".join(charge_values)

    print("To be a valid model, the sum of partial charges must be zero.")
    print("\nCalculating total charge:")
    print(f"{equation_str} = {total_charge:.4f}")

    if abs(total_charge) < 1e-9:
        print("\nThe molecule is electrically neutral. This is a valid set of charges.")
    else:
        print("\nWarning: The molecule is not electrically neutral.")

if __name__ == "__main__":
    propose_methanol_charges()