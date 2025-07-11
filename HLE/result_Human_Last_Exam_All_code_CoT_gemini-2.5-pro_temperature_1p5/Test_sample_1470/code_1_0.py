import sys

def propose_methanol_charges():
    """
    This function proposes a set of reasonable partial charges for methanol,
    verifies its charge neutrality, and prints the results.
    """
    # Based on chemical principles, option D provides the most reasonable charge distribution.
    # C(0.1450), H_methyl(0.0400 x 3), O(-0.6830), H_hydroxyl(0.4180)
    charges = {
        'Carbon': 0.1450,
        'Methyl H1': 0.0400,
        'Methyl H2': 0.0400,
        'Methyl H3': 0.0400,
        'Oxygen': -0.6830,
        'Hydroxyl H': 0.4180
    }

    print("Proposed reasonable partial charges for the atoms of methanol (CH3OH):")
    for atom, charge in charges.items():
        print(f"{atom}:\t{charge:.4f}")

    # Calculate the sum of charges to verify neutrality.
    total_charge = sum(charges.values())

    # Build the equation string as requested
    equation_parts = [f"({c:.4f})" for c in charges.values()]
    equation_str = " + ".join(equation_parts)

    print("\nVerifying that the molecule is charge-neutral:")
    # Using sys.stdout.write to prevent print's default newline and space separator
    sys.stdout.write(f"Sum of charges = {equation_str} = {total_charge:.4f}\n")

    print("\nConclusion: The total charge is zero, making this a valid and chemically reasonable set of parameters.")

if __name__ == '__main__':
    propose_methanol_charges()
